default:
	python3 setup.py install --user

clean:
	rm -fr build
	rm -fr dist
	rm -fr prfect.egg-info
	pip uninstall prfect -y

