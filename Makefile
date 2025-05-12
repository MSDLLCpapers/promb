test: test-full

test-full:
	pytest tests

test-no-deps:
	pytest tests/test_promb_no_deps.py

release:
ifndef VERSION
	$(error "Usage: make release VERSION=0.1.9")
endif
	git checkout main
	git pull
	echo "__version__ = '$(VERSION)'" > promb/__version__.py
	git add promb/__version__.py
	git commit -m "Set version to $(VERSION)"
	git push
	make dist
	twine upload dist/promb-$(VERSION)*

dist:
	python setup.py sdist bdist_wheel	

