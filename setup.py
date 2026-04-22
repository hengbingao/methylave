from setuptools import setup, find_packages

setup(
    name="methylave",
    version="1.0.0",
    description="CG methylation calculator from allc files over BED regions",
    python_requires=">=3.6",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "methylave=methylave.__main__:main",
        ],
    },
)
