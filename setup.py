from setuptools import setup, find_packages

setup(
    name="FastMDAnalysis",
    version="0.1.0",
    packages=find_packages(),
    #packages=find_packages(include=['FastMDAnalysis', 'FastMDAnalysis.*']),
    #package_dir={'FastMDAnalysis': 'FastMDAnalysis'},
    entry_points={
        'console_scripts': [
            'fastmda=FastMDAnalysis.cli:main'
        ]
    },
    install_requires=[
        'mdtraj>=1.9.7',
        'numpy>=1.21.0',
        'matplotlib>=3.5.0',
        'scikit-learn>=1.0.0'
    ],
    author="Adekunle Aina",
    author_email="aaina@csudh.edu",
    description="Fast Molecular Dynamics Trajectory Analysis",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ainaadekunle/FastMDAnalysis",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    include_package_data=True,
)
