from setuptools import setup, find_packages

setup(
    name="famgam",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "tqdm>=4.62.0",
    ],
    author="Jerry Lee",
    author_email="jerrylee9310@gmail.com",
    description="A package for generating family genotypes",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/h2ai/famgam",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
) 