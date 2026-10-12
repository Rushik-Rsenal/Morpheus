import os
from flask import Flask
from . import db

def create_app(test_config=None):
  app = Flask(__name__, instance_relative_config=True)


  app.config.from_mapping(
    SECRET_KEY='dev',
    DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
  )

  # Basically, by defaul test_config = None
  # However, if we pass in something, then we will use the user specified config

  if test_config is None:
    # Load the instance config, if it exists, when not testing
    app.config.from_pyfile("config.py", silent=True)
  else:
    # Load the test config if it exists
    app.config.from_mapping(test_config)

  os.makedirs(app.instance_path, exist_ok=True)

  @app.route('/hello')
  def hello():
    return 'Hello, World!'

  
  db.init_app(app)
  return app