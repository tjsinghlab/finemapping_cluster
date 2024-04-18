import os
import pkg_resources

# define PROJECT_ROOT environment variable
try:
    PROJECT_ROOT = os.path.abspath(
        pkg_resources.resource_filename(pkg_resources.Requirement.parse(__name__), '.')) \
        if 'PROJECT_ROOT' not in os.environ else os.environ['PROJECT_ROOT']
except pkg_resources.DistributionNotFound:
    print('Running without installation. Likely on Spark Cluster.')
    PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

os.environ['PROJECT_ROOT'] = PROJECT_ROOT

# define TRACKER_PATH for use with datatracker
os.environ['TRACKER_PATH'] = os.path.join(
    PROJECT_ROOT, os.path.basename(PROJECT_ROOT), 'db.json') \
    if 'TRACKER_PATH' not in os.environ else os.environ['TRACKER_PATH']

# modify if CLOUD_ROOT is not directly inferred from PROJECT_ROOT
os.environ['CLOUD_ROOT'] = f"gs://{__name__}" \
    if 'CLOUD_ROOT' not in os.environ else os.environ['CLOUD_ROOT']
