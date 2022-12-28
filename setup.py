from pathlib import Path

from skbuild import setup


readme_path = Path("README.md")
long_description = readme_path.read_text(encoding="utf-8")
version_path = Path("VERSION")
version = version_path.read_text(encoding="utf-8").strip()
license_path = Path("LICENSE")
license = license_path.read_text(encoding="utf-8").strip()


setup(  # the cmake build tag used by scikit-build defaults to "Release" 
    name="py-laplace-inversion",
    version=version,
    description="Numerical inversion of Laplace transforms using Den Iseger's algorithm.",
    long_description=long_description,
    url="https://github.com/serban-badila/laplace-inversion",
    author="Serban Badila",
    author_email="badilaserban@gmail.com",
    license=license,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD 3-Clause",
        "Programming Language :: Python :: 3.7",
    ],
    packages=['py_laplace_inversion']
)
