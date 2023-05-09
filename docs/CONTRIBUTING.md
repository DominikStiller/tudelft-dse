# Developer's Guidelines
This file contains guidelines and documentation for developers who want to contribute to this repository. If you only skim these, focus on the **guidelines in bold**.

## Code Guidelines

This project uses the [_Black_ code style](https://github.com/psf/black) for consistency and clarity. Before committing, format the project code by running `black .` in the project directory.  For docstrings and comments, the [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) is used.

Some more guidelines to follow:
* **Use SI units everywhere.**
* **Write a lot of [unit tests](https://docs.python.org/3/library/unittest.html).** This catches errors close to the source and gives you confidence that your code works. If you're using PyCharm, create unit tests for a method by Right-click > Go To > Tests. From the console, all tests can be run using `python -m unittest`.
* Do not use local paths in code that is to be merged into main. This leads to unnecessary changes. Prefer passing paths as command line arguments, which can be set in the PyCharm run configuration. Local paths are okay in experimental code.
* Use [type hints](https://docs.python.org/3/library/typing.html) and [docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) for every method. This helps prevent errors and assists others with using your method properly.
* Ideally, only a single person works on a file at a time to prevent merge conflicts. This requires a certain file structure, avoiding long files and preferring small, specialized files.
* If you start using a package that is not yet in `requirements.txt`, add it with the specific version you are using, so your code works for everyone else as it does for you.



## Git Guidelines
[GitHub Flow](https://docs.github.com/en/get-started/quickstart/github-flow) is used as Git model for this repository. Please become familiar with it before committing. Follow these guidelines:
* **No data should ever be committed to Git.** The repository is for code only. Store any local data files and generated results in the `data/` folder, which is ignored by commits. Data files can be shared through Teams and downloaded by everyone individually.
* **Always start new branch for new feature, do not reuse a branch for multiple features.** A feature should be a single component of about a session's work, if more, split into smaller features. Always update main before starting a new feature branch, then start your feature branch from main.
* **The main branch should always be in a usable state.** This means that pull requests that are not drafts should not contain experimental code. If you want to keep experimental code, but it in `private/`, which is ignored by commits.
* Never commit to the main branch directly (in fact, pushing to main is blocked). Always work on your feature branch, then use a pull request to merge your changes into main.
* If it can be avoided, do not merge feature branches into another. This leads to messy pull requests with much manual labor. Instead, use [cherry-picking](https://gitbetter.substack.com/p/how-to-use-git-cherry-pick-effectively) if you need another branches' code before it has been merged into main.
* Commit often, after finishing a small part of a feature, even if the code does not work fully yet. Since you commit to your own feature branch, no one else is affected, and you can keep a history of your changes in case something goes wrong.
* Use descriptive commit messages (not just "fix bugs"). Follow [these guidelines](https://gist.github.com/robertpainsi/b632364184e70900af4ab688decf6f53). Use the imperative ("add" instead of "added") for verbs.
* Use pull requests as platform for discussions and questions. You can open a pull request even if your code is not done yet. Tagging people in pull requests to get feedback on work-in-progress code is explicitly encouraged. Once you're done, the pull request will be approved for merging into main.
