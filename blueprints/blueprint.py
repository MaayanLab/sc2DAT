from flask import Blueprint, url_for
import os
import s3fs
from dotenv import load_dotenv
load_dotenv()
blueprint = Blueprint('blueprint', __name__)


@blueprint.route('numruns')
def numruns():
  key = os.getenv('APPYTER_DATA_DIR').split('#?key=')[-1].split('&')[0]
  secret = os.getenv('APPYTER_DATA_DIR').split('&secret=')[-1]
  s3 = s3fs.S3FileSystem(key = key, secret = secret)
  print(len(s3.ls('s3://multiomics2paper/output')))
  n_runs = len(s3.ls('s3://multiomics2paper/output'))
  return {'num_runs': n_runs}