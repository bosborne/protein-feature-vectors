from setuptools import setup, find_packages

setup(
    name="protein-feature-vectors",
    version="1.0",
    keywords="protein-feature-vectors",
    description="Code to produce fixed-length feature vectors from protein sequences based on iFeatureOmega",
    license="MIT License",
    url="https://github.com/bosborne/protein-feature-vectors",
    author="Brian Osborne",
    author_email="bosborne@bioteam.net",
    packages=find_packages("src"),
    package_dir={"": "src"},
    package_data={
        "": ["*.txt"],
        "protein-feature-vectors": [
            "data/*.txt",
            "data/*.data",
            "data_examples/*.csv",
            "data_examples/*.fa",
            "*.json",
        ],
    },
    platforms="any",
    install_requires=[
        "numpy>=1.21.4",
        "pandas>=1.3.4",
        "scipy>=1.7.3",
        "biopython>=1.6",
    ],
)
