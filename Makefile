default:
	python3 setup.py install --user

clean:
	rm -fr build
	rm -fr dist
	rm -fr slippery.egg-info
	pip uninstall slipery -y

