from flask import Blueprint, make_response
import os
import s3fs
from dotenv import load_dotenv
load_dotenv()
blueprint = Blueprint('blueprint', __name__)

@blueprint.route('numruns')
def numruns():
  if 'APPYTER_DATA_DIR' not in os.environ or '#?key=' not in os.getenv('APPYTER_DATA_DIR') or '&secret=' not in os.getenv('APPYTER_DATA_DIR'):
    return {'num_runs': 0}
  key = os.getenv('APPYTER_DATA_DIR').split('#?key=')[-1].split('&')[0]
  secret = os.getenv('APPYTER_DATA_DIR').split('&secret=')[-1]
  s3 = s3fs.S3FileSystem(key = key, secret = secret)
  n_runs = len(s3.ls('s3://ct2targets/output'))
  response = make_response({'num_runs': n_runs})
  response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
  response.headers['Pragma'] = 'no-cache'
  response.headers['Expires'] = '0'
  return response