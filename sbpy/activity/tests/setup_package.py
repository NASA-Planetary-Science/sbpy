import os


def get_package_data():
    return {'sbpy.activity.tests': [os.path.join('data', '*.fits')]}
