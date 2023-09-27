from setuptools import setup, find_packages
import os

setup(
    name='eboss_qso',
    version='0.0.1',
    author='Nick Hand',
    package_data = {'eboss_qso': ['fits/nersc/job.template.sh']},
    packages=find_packages(),
    description=("Python code to analyze primordial non-Gaussianity "
                 "from the eBOSS QSO sample"),
    install_requires=['numpy', 'scipy'],
    entry_points={
        'console_scripts': [
            'eboss-echo-hash = eboss_qso.measurements.utils:echo_hash',
            'eboss-qso-fit = eboss_qso.fits.driver:__main__',
            'eboss-genreports = eboss_qso.fits.reports:_generate_all',
            'eboss-gentoc = eboss_qso.fits.reports:generate_report_toc',
            'eboss-sync = eboss_qso.sync:sync_eboss',
            'eboss-fix-paths = eboss_qso.fits.results:fix_result_paths',
            ]
    }
)
