from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename
from flaskr.auth import login_required
from flaskr.db import get_db
import os
import sys
import subprocess
import tempfile
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import pyplot as plt
import sys
import numpy as np
import io
import base64

ALLOWED_DATA_EXTENSIONS = {'csv', 'json', 'txt', 'xlsx'}
ALLOWED_SCRIPT_EXTENSIONS = {'py'}

def allowed_file(filename, allowed_set):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in allowed_set


bp = Blueprint('blog', __name__)

plot_url = None

@bp.route('/')
def index():
    global plot_url

    db = get_db()
    posts = db.execute(
        'SELECT p.id, title, body, created, author_id, username'
        ' FROM post p JOIN user u ON p.author_id = u.id'
        ' ORDER BY created DESC'
    ).fetchall()

    return render_template('blog/index.html', posts=posts, plot_url=plot_url)

@bp.route('/create', methods=('GET', 'POST'))
@login_required
def create():
    global plot_url

    if request.method == 'POST':
        title = request.form['title']
        body = request.form['body']

        # Check if the file is there
        if 'script' not in request.files:
            flash('No file part')
            return redirect(request.url)

        script = request.files['script']
        dataset = request.files['dataset']


        if script.filename == '':
            flash('No selected script')
            return redirect(request.url)
        elif dataset.filename == '':
            flash("No selected dataset")
            return redirect(request.url)
        if script and allowed_file(script.filename, ALLOWED_SCRIPT_EXTENSIONS) and allowed_file(dataset.filename, ALLOWED_DATA_EXTENSIONS):

            # We store the dataset
            db = get_db()

            # Check the file extension
            if dataset.filename.rsplit('.', 1)[1].lower() == 'xlsx':
                # Read the excel file
                df = pd.read_excel(dataset, sheet_name=0)

                # Update the params for better looking graphs
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


                # Set the index

                # Look at the feature
                stressor = request.form["stressor-name"]
                feature = request.form["feature-name"]



                # Develop the image
                layout = """
                    000000000000
                    111122223333
                    444455556666
                """

                fig, axes = plt.subplot_mosaic(layout, figsize=(30, 20), gridspec_kw={"wspace": 1, "hspace": 0.5})

                # Add text title to the top left (Area '0')
                axes['0'].axis('off') # Hide grid lines
                axes['0'].text(0.5, 0.5, f"{feature} ({stressor})", 
                            fontsize=50, ha='center', va='center')

                # Log-fold change function
                def well_level_data_analysis_log_fold_change(df, col_name, timepoint, stressor, letter, axes):
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

                    for batch in batches:
                        selected_batch = df_extracted.loc[batch]

                        # Get the control rows
                        control_rows = selected_batch[selected_batch["Celltype"] == "Control"]

                        # Get the average
                        avg_df = control_rows.groupby(["Celltype"])[col_name].mean()

                        # Store this in the dictionary created
                        control_avg_dict[batch] = avg_df["Control"]

                        # Divide each row in the selected_batch by the mean of controls
                        selected_batch_normalised = selected_batch
                        selected_batch_normalised[col_name] = np.log(selected_batch[col_name]/control_avg_dict[batch])


                        # # Now get the average for cellines
                        cellline_df = selected_batch_normalised.groupby(["Cellline ID", "Celltype"])[col_name].mean()

                        # # Append to the list
                        cellline_df_list.append(cellline_df)

                    # Combine the cellline_df_list
                    combined_cellline_df = pd.concat(cellline_df_list, axis=0)
                    combined_cellline_df = combined_cellline_df.reset_index()
                    combined_cellline_df = combined_cellline_df.sort_values(by="Celltype")

                    # Visualise the stuff 
                    combined_cellline_df_excluded = combined_cellline_df[combined_cellline_df["Celltype"] != "Control"]
                    sns.barplot(data=combined_cellline_df_excluded, y="Cellline ID", x=col_name, hue="Celltype", ax=axes[str(timepoint)]).set_title(f"{letter}) log-fold change of {col_name} {timepoint}")
                    plt.figure(figsize=(8, 10),)

                # Plot your graphs
                def return_img():
                    well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 1, "DMSO", "a", axes)
                    well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 2, "DMSO", "b", axes)
                    well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 3, "DMSO", "c", axes)
                    well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 4, "DMSO", "d", axes)
                    well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 5, "DMSO", "e", axes)
                    well_level_data_analysis_log_fold_change(df, "mito_TMRM intensity within MTG", 6, "DMSO", "f", axes)
                    return fig

                # Return the image
                fig = return_img()

                # Save the figure to a BytesIO buffer
                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format="png", bbox_inches="tight")
                img_buffer.seek(0)

                # Encode the buffer content to base64
                plot_url = base64.b64encode(img_buffer.getvalue()).decode("utf8")

                # Clean up memory by closing the figure
                plt.close(fig)

                

                error = None

                if not title:
                    error = 'Title is required.'

                if error is not None:
                    flash(error)
                else:
                    db = get_db()
                    db.execute(
                        'INSERT INTO post (title, body, author_id, script_name, script_body)'
                        ' VALUES (?, ?, ?, ?, ?)',
                        (title, body, g.user['id'], script.filename, script.read())
                    )
                    db.commit()
                    # 5. Pass the base64 string to the HTML template
                    return redirect(url_for('blog.index'))
                    # return render_template("blog/index.html", plot_url=plot_url)
            else:
                print("Bro wrong datatype")



    return render_template('blog/create.html')

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






@bp.route('/<int:id>/use', methods=('GET', 'POST'))
@login_required
def use(id):
    output = None
    error = None

    if request.method == 'POST':
        dataset = request.files.get('dataset')
        script = request.files.get('script')
        print(script)

        if not dataset or not script or dataset.filename == '' or script.filename == '':
            error = "Both a dataset and a Python script are required."
        elif not allowed_file(dataset.filename, ALLOWED_DATA_EXTENSIONS):
            error = f"Invalid dataset format. Allowed: {', '.join(ALLOWED_DATA_EXTENSIONS)}"
        elif not allowed_file(script.filename, ALLOWED_SCRIPT_EXTENSIONS):
            error = "Invalid script format. Only .py files are allowed."
        else:
            # Create an isolated temporary directory that auto-deletes when done
            with tempfile.TemporaryDirectory() as temp_dir:
                data_filename = secure_filename(dataset.filename)
                script_filename = secure_filename(script.filename)

                data_path = os.path.join(temp_dir, data_filename)
                script_path = os.path.join(temp_dir, script_filename)

                dataset.save(data_path)
                script.save(script_path)

                try:
                    # Execute the uploaded script using the current Python interpreter
                    # Passes the dataset file path as the first argument (sys.argv[1])
                    result = subprocess.run(
                        [sys.executable, script_path, data_path],
                        capture_output=True,
                        text=True,
                        timeout=30,  # Prevent infinite loops
                        cwd=temp_dir # Run inside the temp directory
                    )
                    output = result.stdout
                    if result.stderr:
                        error = result.stderr
                except subprocess.TimeoutExpired:
                    error = "Execution timed out (30-second limit exceeded)."
                except Exception as e:
                    error = f"Failed to execute script: {str(e)}"

    return redirect(url_for('blog.index'))