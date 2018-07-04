python3 setup.py sdist
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
sudo pip3 uninstall kmalgorithm
sudo pip3 install --index-url https://test.pypi.org/simple/ kmalgorithm

