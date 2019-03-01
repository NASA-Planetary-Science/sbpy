import os


def get_package_data():
    return {'sbpy.photometry.tests': [os.path.join('data', '*.fits'),
                                      os.path.join('data', '*.txt')]}
