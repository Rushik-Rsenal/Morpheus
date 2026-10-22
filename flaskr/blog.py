from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename
from flaskr.auth import login_required
from flaskr.db import get_db
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import pyplot as plt
import numpy as np
import io
import base64
import urllib.parse

ALLOWED_DATA_EXTENSIONS = {'csv', 'json', 'txt', 'xlsx'}
ALLOWED_SCRIPT_EXTENSIONS = {'py'}

def allowed_file(filename, allowed_set):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in allowed_set


bp = Blueprint('blog', __name__)

initial_dataframe = None

@bp.route('/')
def index():
    # Get the db
    db = get_db()

    # Get the posts
    posts = db.execute(
        'SELECT p.id, title, body, created, author_id, username'
        ' FROM post p JOIN user u ON p.author_id = u.id'
        ' ORDER BY created DESC'
    ).fetchall()

    # Get the images
    images = db.execute(
        'select plot_url, post_id'
        ' from image_table i join post on post.id = i.post_id'
    ).fetchall()

    return render_template('blog/index.html', posts=posts, images=images)

# Create a workflow
@bp.route('/create', methods=('GET', 'POST'))
@login_required
def create():
    # This is a global variable used to store the user's uploaded file
    global initial_dataframe

    # Look at the POST request
    if request.method == 'POST':
        # Get the title and the body
        title = request.form['title']
        body = request.form['body']

        # Get the dataset
        dataset = request.files['dataset']

        # Get the dataset filename
        if dataset.filename == '':
            flash("No selected dataset")
            return redirect(url_for("blog.create"))

        # Get the allowed data extensions
        if allowed_file(dataset.filename, ALLOWED_DATA_EXTENSIONS):
            # Check the file extension
            if dataset.filename.rsplit('.', 1)[1].lower() == 'xlsx':
                # Read the excel file
                df = pd.read_excel(dataset, sheet_name=0)
                initial_dataframe = df

                # Get only the features and not other well information 
                features = initial_dataframe.drop(
                    columns=["Celltype", 
                             "Cellline ID",
                              "Batch",
                              "Stressor",
                              "Row",
                              "Column",
                              "Timepoint",
                              "cell_number"]
                ).columns
                columns = features

                # Get the timepoints 
                timepoints = initial_dataframe["Timepoint"].unique() 

                # Get the cell types
                cell_types = initial_dataframe["Celltype"].unique()

                # Get the cellline IDs
                cell_lines = initial_dataframe["Cellline ID"].unique()

                # Get the stressors
                stressors = initial_dataframe["Stressor"].unique()

                # Get the batches
                batches = initial_dataframe["Batch"].unique()

                # Basic Error handling
                error = None

                if not title:
                    error = 'Title is required.'
                
                if error is not None:
                    flash(error)
                else:
                    # Insert the post related stuff 
                    db = get_db()
                    db.execute(
                        'INSERT INTO post (title, body, author_id, script_name, script_body)'
                        ' VALUES (?, ?, ?, ?, ?)',
                        (title, body, g.user['id'], "script.filename", "script.read()")
                    )
                    db.commit()

                    # Redirect to the next screen
                    query_string = urllib.parse.urlencode( {'timepoints': timepoints, 
                                                           "cell_types": cell_types,
                                                           "stressors": stressors,
                                                           "batches": batches,
                                                           "columns": columns,
                                                           "cell_lines": cell_lines}, doseq=True)
    
                    # Construct the full destination URL and redirect
                    target_url = f"{url_for('blog.next_screen')}?{query_string}"
                    return redirect(target_url)
            else:
                print("Bro wrong datatype")
    return render_template('blog/create.html')

# Log-fold change function
def well_level_data_analysis_log_fold_change(df, col_name, timepoint, stressor, letter, axes):
    #  Update the params for better looking graphs
    plt.rcParams.update(
                        {
                            "font.size": 11,
                            "font.family": 'serif',
                            "font.serif": ['Liberation Serif'],
                            "mathtext.fontset": 'stix',
                            "axes.labelsize": 12,
                            "xtick.labelsize": 10,
                            "ytick.labelsize": 10,
                            "legend.fontsize": 10,
                            "figure.figsize": (6,4),
                            "figure.dpi": 300,
                            "savefig.dpi": 300,
                            "savefig.bbox": "tight",
                            "axes.linewidth": 1.2,
                            "grid.linewidth": 0.8,
                            "lines.linewidth": 2,
                        }
    )

    # Extract rows from a particular timepoint
    df_init = df
    df_extracted = df_init[df_init["Timepoint"] == timepoint]
    df_extracted = df_extracted[df_extracted["Stressor"] == stressor]

    # Extract the required columns
    df_extracted = df_extracted[["Celltype", "Cellline ID","Batch", col_name]]

    # Calculate the average of Control for each Batch
    batches = df_extracted["Batch"].unique()
    df_extracted = df_extracted.set_index("Batch")

    # Control_AVG
    control_avg_dict = {}
    cellline_df_list = []

    # For each batch, calculate the cell-line average
    for batch in batches:
        # Get the batch
        selected_batch = df_extracted.loc[batch]

        # Error check
        if len(selected_batch) == 0:
            flash("This messed up!")

        # Get the control rows
        control_rows = selected_batch[selected_batch["Celltype"] == "Control"]

        # Get the average
        avg_df = control_rows.groupby(["Celltype"])[col_name].mean()

        # Store this in the dictionary created
        control_avg_dict[batch] = avg_df["Control"]

        # Divide each row in the selected_batch by the mean of controls
        selected_batch_normalised = selected_batch
        selected_batch_normalised[col_name] = np.log(selected_batch[col_name]/control_avg_dict[batch])


        # Now get the average for cellines
        cellline_df = selected_batch_normalised.groupby(["Cellline ID", "Celltype"])[col_name].mean()
        cellline_df = cellline_df.reset_index()

        # Append to the list
        cellline_df_list.append(cellline_df)

    if len(cellline_df_list) == 0:
        flash("Bro something messed up!")

    # Combine the cellline_df_list
    combined_cellline_df = pd.concat(cellline_df_list, axis=0)
    combined_cellline_df = combined_cellline_df.reset_index()
    combined_cellline_df = combined_cellline_df.sort_values(by="Celltype")

    # Visualise the stuff 
    combined_cellline_df_excluded = combined_cellline_df[combined_cellline_df["Celltype"] != "Control"]
    sns.barplot(data=combined_cellline_df_excluded, y="Cellline ID", x=col_name, hue="Celltype", ax=axes[str(timepoint)]).set_title(f"{letter}) log-fold change of {col_name} {timepoint}")
    plt.figure(figsize=(8, 10),)



@bp.route('/next_screen', methods=('GET', 'POST'))
@login_required
def next_screen():
    global images, initial_dataframe
    timepoints = request.args.getlist("timepoints")
    cell_types = request.args.getlist("cell_types")
    stressors = request.args.getlist("stressors")
    columns = request.args.getlist("columns")
    batches = request.args.getlist("batches")
    cell_lines = request.args.getlist("cell_lines")

    # POST Method check
    if request.method == "POST":
        # Get the features
        features = initial_dataframe.drop(columns=["Celltype", "Cellline ID", "Batch", "Stressor", "Row", "Column", "Timepoint", "cell_number"]).columns
        
        # Get the timepoints 
        timepoints = initial_dataframe["Timepoint"].unique() 
        print("Timepoints: ", timepoints)
        
        # Get the cell types
        cell_types = initial_dataframe["Celltype"].unique()
        print("Cell Types: ", cell_types)

         # Get the cell lines
        cell_lines = initial_dataframe["Cellline ID"].unique()
        print("Cell Lines: ", cell_lines)

        # Get the stressors
        stressors = initial_dataframe["Stressor"].unique()
        print("Stressors: ", stressors)

        # Get the batches
        batches = initial_dataframe["Batch"].unique()
        print("Batches: ", batches)

        # Get the checked features
        selected_features = []

        # Go over the form
        for feature in features:
            if request.form.get(feature):
                selected_features.append(feature)

        print("Features: ", selected_features)

        # Get the checked timepoints
        selected_timepoints = []

        # Go over the form
        for timepoint in timepoints:
            if request.form.get(f"T{timepoint}"):
                selected_timepoints.append(timepoint)

        print("Timepoints: ", selected_timepoints)

        # Get the checked stressor
        selected_stressors = []

        # Go over the form
        for stressor in stressors:
            if request.form.get(stressor):
                selected_stressors.append(stressor)

        print("stressors: ", selected_stressors)

        # Get the checked batches
        selected_batches = []
        
        # Go over the form
        for batch in batches:
            if request.form.get(batch):
                selected_batches.append(batch)

        print("batches: ", selected_batches)

        # Get the checked cell_types
        selected_cell_types = []
        # Go over the form
        for cell_type in cell_types:
            if request.form.get(cell_type):
                selected_cell_types.append(cell_type)

        print("cell_types ", selected_cell_types)

        # Get the checked cell lines
        selected_cell_lines = []

        # Go over the form
        for cell_line in cell_lines:
            if request.form.get(cell_line):
                selected_cell_lines.append(cell_line)

        print("cell_lines", selected_cell_lines)

        # Select the necessary columns
        selected_columns = selected_features


        arr = ["Timepoint", "Batch", "Cellline ID", "Celltype", "Stressor"]
        for i in arr:
            selected_columns.append(i)
        
        print("selected columns ", selected_columns)
        selected_df = initial_dataframe[selected_columns]

        print("DF before")
        print(selected_df.head(2))

        # Now do great things
        # Filter on the necessary timepoints
        selected_df = selected_df[selected_df["Timepoint"].isin(selected_timepoints)]
        print("After timepoint")
        print(selected_df.head(2))


        # Filter on the necessary celltypes
        selected_df = selected_df[selected_df["Celltype"].isin(selected_cell_types)]
        print("After cell type")
        print(selected_df.head(2))

        # Filter on the necessary celllines
        selected_df = selected_df[selected_df["Cellline ID"].isin(selected_cell_lines)]

        # Filter on the necessary batches
        selected_df = selected_df[selected_df["Batch"].isin(selected_batches)]
        print("After batch")
        print(selected_df.head(2))

        # For each stressor
            # For each feature
                # Produce a nice dashboard
        images = []
        for stressor in selected_stressors:
            for feature in selected_features:
                # Develop the image
                layout = """
                        00000
                        11111
                        11111
                        """

                fig, axes = plt.subplot_mosaic(layout, figsize=(30, 20), gridspec_kw={"wspace": 1, "hspace": 0.5})

                # Title
                axes['0'].axis('off') 
                axes['0'].text(0.5, 0.5, f"{feature} ({stressor})", fontsize=50, ha='center', va='center')

            
                # Plot your graphs
                count = 0
                for timepoint in selected_timepoints:
                    well_level_data_analysis_log_fold_change(selected_df,feature , timepoint, stressor, f"{count}", axes)
                    count += 1
                
                # Save the figure to a BytesIO buffer
                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format="png", bbox_inches="tight")
                img_buffer.seek(0)

                # Encode the buffer content to base64
                plot_url = base64.b64encode(img_buffer.getvalue()).decode("utf8")

                # Clean up memory by closing the figure
                db = get_db()

                # First get the required post id
                post_id = db.execute(
                    'SELECT id FROM post WHERE author_id = ? ORDER BY created DESC LIMIT 1',
                    (g.user['id'],)
                ).fetchone()['id']
                 
                # Then, use the that post_id and add it to the images table
                db.execute(
                    'insert into image_table(post_id, plot_url)'
                    ' values (?, ?)'
                    ,
                    (post_id, plot_url)
                )
                db.commit()
                images.append(plot_url)

                plt.close(fig)
                db = get_db()
                posts = db.execute(
                        'SELECT p.id, title, body, created, author_id, username'
                        ' FROM post p JOIN user u ON p.author_id = u.id'
                        ' ORDER BY created DESC'
                ).fetchall()
                
                images = db.execute(
                        'select plot_url, post_id'
                        ' from image_table i join post on post.id = i.post_id'
                ).fetchall()
                
                return redirect(url_for('blog.index'))
    
    
    return render_template('blog/next_screen.html', batches=batches, columns=columns, stressors=stressors, cell_types=cell_types, timepoints=timepoints, cell_lines=cell_lines)


# Box plot
# Create a workflow
@bp.route('/create_box_plot', methods=('GET', 'POST'))
@login_required
def create_box_plot():
    global initial_dataframe

    if request.method == 'POST':
        # Get the title and the body
        title = request.form['title']
        body = request.form['body']

        # Get the dataset
        dataset = request.files['dataset']

        # Get the dataset filename
        if dataset.filename == '':
            flash("No selected dataset")
            return redirect(url_for("blog.create_box_plot"))

        # Get the allowed data extensions
        if allowed_file(dataset.filename, ALLOWED_DATA_EXTENSIONS):
            # Check the file extension
            if dataset.filename.rsplit('.', 1)[1].lower() == 'xlsx':
                # Read the excel file
                df = pd.read_excel(dataset, sheet_name=0)
                initial_dataframe = df

                # Get only the features and not other well information 
                features = initial_dataframe.drop(
                    columns=["Celltype", 
                             "Cellline ID",
                              "Batch",
                              "Stressor",
                              "Row",
                              "Column",
                              "Timepoint",
                              "cell_number"]
                ).columns
                columns = features

                # Get the timepoints 
                timepoints = initial_dataframe["Timepoint"].unique() 

                # Get the cell types
                cell_types = initial_dataframe["Celltype"].unique()

                # Get the stressors
                stressors = initial_dataframe["Stressor"].unique()

                # Get the batches
                batches = initial_dataframe["Batch"].unique()

                # Get the cell lines
                cell_lines = initial_dataframe["Cellline ID"].unique()

                # Basic Error handling
                error = None

                if not title:
                    error = 'Title is required.'
                
                if error is not None:
                    flash(error)
                else:
                    # Insert the post related stuff 
                    db = get_db()
                    db.execute(
                        'INSERT INTO post (title, body, author_id, script_name, script_body)'
                        ' VALUES (?, ?, ?, ?, ?)',
                        (title, body, g.user['id'], "script.filename", "script.read()")
                    )
                    db.commit()

                    # Redirect to the next screen
                    query_string = urllib.parse.urlencode({'timepoints': timepoints, 
                                                           "cell_types": cell_types,
                                                           "cell_lines": cell_lines,
                                                           "stressors": stressors,
                                                           "batches": batches,
                                                           "columns": columns}, doseq=True)
                    
                    # Construct the full destination URL and redirect
                    target_url = f"{url_for('blog.next_screen_box_plot')}?{query_string}"
                    return redirect(target_url)
            else:
                print("Bro wrong datatype")
    return render_template('blog/create.html')

# Box plot
def well_level_box_plot(df, col_name, timepoint, stressor, letter, axes):
    #  Update the params for better looking graphs
    plt.rcParams.update(
                        {
                            "font.size": 13,
                            "font.family": 'serif',
                            "font.serif": ['Liberation Serif'],
                            "mathtext.fontset": 'stix',
                            "axes.labelsize": 12,
                            "xtick.labelsize": 10,
                            "ytick.labelsize": 10,
                            "legend.fontsize": 10,
                            "figure.figsize": (6,4),
                            "figure.dpi": 300,
                            "savefig.dpi": 300,
                            "savefig.bbox": "tight",
                            "axes.linewidth": 1.2,
                            "grid.linewidth": 0.8,
                            "lines.linewidth": 2,
                        }
    )

    # Extract rows from a particular timepoint
    df_init = df
    df_extracted = df_init[df_init["Timepoint"] == timepoint]
    df_extracted = df_extracted[df_extracted["Stressor"] == stressor]

    # Extract the required columns
    df_extracted = df_extracted[["Celltype", "Cellline ID","Batch", col_name]]

    # Calculate the average of Control for each Batch
    batches = df_extracted["Batch"].unique()
    df_extracted = df_extracted.set_index("Batch")

    # Control_AVG
    control_avg_dict = {}
    cellline_df_list = []

    # For each batch, calculate the cell-line average
    for batch in batches:
        # Get the batch
        selected_batch = df_extracted.loc[batch]

        # Error check
        if len(selected_batch) == 0:
            flash("This messed up!")

        # Get the control rows
        control_rows = selected_batch[selected_batch["Celltype"] == "Control"]

        # Get the average
        avg_df = control_rows.groupby(["Celltype"])[col_name].mean()

        # Store this in the dictionary created
        control_avg_dict[batch] = avg_df["Control"]

        # Divide each row in the selected_batch by the mean of controls
        selected_batch_normalised = selected_batch
        selected_batch_normalised[col_name] = selected_batch[col_name]/control_avg_dict[batch]


        # Now get the average for cellines
        cellline_df = selected_batch_normalised.groupby(["Cellline ID", "Celltype"])[col_name].mean()
        cellline_df = cellline_df.reset_index()

        # Append to the list
        cellline_df_list.append(cellline_df)

    if len(cellline_df_list) == 0:
        flash("Bro something messed up!")

    # Combine the cellline_df_list
    combined_cellline_df = pd.concat(cellline_df_list, axis=0)
    combined_cellline_df = combined_cellline_df.reset_index()
    combined_cellline_df = combined_cellline_df.sort_values(by="Celltype")

    # Visualise the stuff 
    combined_cellline_df_excluded = combined_cellline_df
    sns.boxplot(data=combined_cellline_df_excluded, 
                y="Celltype",
                x=col_name,
                hue="Celltype",
                ax=axes['1']).set_title(f"{letter}) {col_name} at {timepoint}")
    plt.figure(figsize=(8, 10),)



@bp.route('/next_screen_box_plot', methods=('GET', 'POST'))
@login_required
def next_screen_box_plot():
    global images, initial_dataframe
    timepoints = request.args.getlist("timepoints")
    cell_types = request.args.getlist("cell_types")
    stressors = request.args.getlist("stressors")
    columns = request.args.getlist("columns")
    batches = request.args.getlist("batches")
    cell_lines = request.args.getlist("cell_lines")

    # POST Method check
    if request.method == "POST":
        # Get the features
        features = initial_dataframe.drop(columns=["Celltype", "Cellline ID", "Batch", "Stressor", "Row", "Column", "Timepoint", "cell_number"]).columns
        
        # Get the timepoints 
        timepoints = initial_dataframe["Timepoint"].unique() 
        print("Timepoints: ", timepoints)
        
        # Get the cell types
        cell_types = initial_dataframe["Celltype"].unique()
        print("Cell Types: ", cell_types)

        # Get the stressors
        stressors = initial_dataframe["Stressor"].unique()
        print("Stressors: ", stressors)

        # Get the batches
        batches = initial_dataframe["Batch"].unique()
        print("Batches: ", batches)

        # Get the batches
        cell_lines = initial_dataframe["Cellline ID"].unique()
        print("cell_lines: ", cell_lines)



        # Get the checked features
        selected_features = []

        # Go over the form
        for feature in features:
            if request.form.get(feature):
                selected_features.append(feature)

        print("Features: ", selected_features)

        # Get the checked timepoints
        selected_timepoints = []

        # Go over the form
        for timepoint in timepoints:
            if request.form.get(f"T{timepoint}"):
                selected_timepoints.append(timepoint)

        print("Timepoints: ", selected_timepoints)

        # Get the checked stressor
        selected_stressors = []

        # Go over the form
        for stressor in stressors:
            if request.form.get(stressor):
                selected_stressors.append(stressor)

        print("stressors: ", selected_stressors)

        # Get the checked batches
        selected_batches = []
        
        # Go over the form
        for batch in batches:
            if request.form.get(batch):
                selected_batches.append(batch)

        print("batches: ", selected_batches)

        # Get the checked cell_types
        selected_cell_types = []
        # Go over the form
        for cell_type in cell_types:
            if request.form.get(cell_type):
                selected_cell_types.append(cell_type)

        print("cell_types ", selected_cell_types)


        selected_cell_lines = []

        # Go over the form
        for cell_line in cell_lines:
            if request.form.get(cell_line):
                selected_cell_lines.append(cell_line)

        print("cell_lines", selected_cell_lines)

        # Select the necessary columns
        selected_columns = selected_features


        arr = ["Timepoint", "Batch", "Cellline ID", "Celltype", "Stressor"]
        for i in arr:
            selected_columns.append(i)
        
        print("selected columns ", selected_columns)
        selected_df = initial_dataframe[selected_columns]

        print("DF before")
        print(selected_df.head(2))

        # Now do great things
        # Filter on the necessary timepoints
        selected_df = selected_df[selected_df["Timepoint"].isin(selected_timepoints)]
        print("After timepoint")
        print(selected_df.head(2))


        # Filter on the necessary celltypes
        selected_df = selected_df[selected_df["Celltype"].isin(selected_cell_types)]
        print("After cell type")
        print(selected_df.head(2))

        # Filter on the necessary celltypes
        selected_df = selected_df[selected_df["Cellline ID"].isin(selected_cell_lines)]
        print("After cellline filtering")
        print(selected_df.head(2))

        # Filter on the necessary batches
        selected_df = selected_df[selected_df["Batch"].isin(selected_batches)]
        print("After batch")
        print(selected_df.head(2))

        # For each stressor
            # For each feature
                # Produce a nice dashboard
        images = []
        for stressor in selected_stressors:
            for feature in selected_features:
                # Develop the image
                layout = """
                        00000
                        11111
                        11111
                        """

                fig, axes = plt.subplot_mosaic(layout, figsize=(30, 20), gridspec_kw={"wspace": 1, "hspace": 0.5})

                # Title
                axes['0'].axis('off') 
                axes['0'].text(0.5, 0.5, f"{feature} ({stressor})", fontsize=50, ha='center', va='center')

            
                # Plot your graphs
                count = 0
                for timepoint in selected_timepoints:
                    well_level_box_plot(selected_df,feature , timepoint, stressor, f"{count}", axes)
                    count += 1
                
                # Save the figure to a BytesIO buffer
                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format="png", bbox_inches="tight")
                img_buffer.seek(0)

                # Encode the buffer content to base64
                plot_url = base64.b64encode(img_buffer.getvalue()).decode("utf8")

                # Clean up memory by closing the figure
                db = get_db()

                # First get the required post id
                post_id = db.execute(
                    'SELECT id FROM post WHERE author_id = ? ORDER BY created DESC LIMIT 1',
                    (g.user['id'],)
                ).fetchone()['id']
                 
                # Then, use the that post_id and add it to the images table
                db.execute(
                    'insert into image_table(post_id, plot_url)'
                    ' values (?, ?)'
                    ,
                    (post_id, plot_url)
                )
                db.commit()
                images.append(plot_url)

                plt.close(fig)
                db = get_db()
                posts = db.execute(
                        'SELECT p.id, title, body, created, author_id, username'
                        ' FROM post p JOIN user u ON p.author_id = u.id'
                        ' ORDER BY created DESC'
                ).fetchall()
                
                images = db.execute(
                        'select plot_url, post_id'
                        ' from image_table i join post on post.id = i.post_id'
                ).fetchall()
                
                # return render_template('blog/index.html', posts=posts, images=images)
                return redirect(url_for('blog.index'))
    
    
    return render_template('blog/next_screen.html', batches=batches, columns=columns, stressors=stressors, cell_types=cell_types, timepoints=timepoints, cell_lines=cell_lines)


def get_post(id, check_author=True):
    post = get_db().execute(
        'SELECT p.id, title, body, created, author_id, username, script'
        ' FROM post p JOIN user u ON p.author_id = u.id'
        ' WHERE p.id = ?',
        (id,)
    ).fetchone()

    if post is None:
        abort(404, f"Post id {id} doesn't exist.")

    if check_author and post['author_id'] != g.user['id']:
        abort(403)

    return post

@bp.route('/<int:id>/update', methods=('GET', 'POST'))
@login_required
def update(id):
    post = get_post(id)

    if request.method == 'POST':
        title = request.form['title']
        body = request.form['body']
        script = post["script"]
        error = None

        if not title:
            error = 'Title is required.'

        if error is not None:
            flash(error)
        else:
            db = get_db()
            db.execute(
                'UPDATE post SET title = ?, body = ?, script = ? WHERE id = ?',
                (title, body, str(script), id)
            )
            db.commit()
            return redirect(url_for('blog.index'))

    return render_template('blog/update.html', post=post)


@bp.route('/<int:id>/delete', methods=('POST',))
@login_required
def delete(id):
    get_post(id)
    db = get_db()
    db.execute('DELETE FROM post WHERE id = ?', (id,))
    db.commit()
    return redirect(url_for('blog.index'))