# CONTRIBUTING

By contributing to SpatialPy, you agree to the [Code of Conduct](/CODE_OF_CONDUCT.md) and that the contents of this project are distributed under the [GNU GPL v3](/LICENSE). Please review both of these documents. If after reading you are interested in contributing, please follow these guidelines:  

### 1. Making Changes:  
a) Create a fork from this repository ('https://github.com/spatialpy/SpatialPy')  
b) Create a new branch.  
c) Please use a meaningful name, and include the issue # if applicable.  
d) For bug fixes, branch from 'master'  
e) For new features, branch from 'develop'  
f) Be sure to document your code.  
g) Eliminate any commented-out or dead code.  
h) If you are creating a new file, prepend the [COPYRIGHT](/COPYRIGHT) notice to the beginning of the file.  
  
### 2. Submitting a Pull Request:  
a) If changes are bug/hotfix, make a pull request to 'master'.  
b) If other changes, or new features are being added, make a pull request to 'develop'.  
c) Include a list of changes in the PR description.  
d) Provide a usage guide/how-to for new features in the PR description.  
  
***Do NOT merge your own pull request.  Once the request is submitted, the code must be reviewed and merged by someone else.***  
  
### 3. Merging a Pull Request:
a) Verify correct merge destination ('master' for hotfix/bugs, 'develop' for features/changes).  
b) Review code for logic, consistency, documentation, and commented-out or dead sections.
c) Verify that unit tests are provided for the new code, and that they accurately test the new feature/changes.  
d) Make sure that all unit tests in Travis CI pass before merging the changes.  
e) Merge!  
