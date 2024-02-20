# docs/Makefile.am

doc_DATA  = audiowmark.pdf audiowmark.html

GRAPHVIZ_PY = graphviz.py

audiowmark.pdf: audiowmark.md $(GRAPHVIZ_PY)
	pandoc -F $(GRAPHVIZ_PY) -V papersize:a4 -V geometry:margin=2cm $< -o $@

audiowmark.html: audiowmark.md $(GRAPHVIZ_PY)
	pandoc -F $(GRAPHVIZ_PY) $< -o $@

# DEPS: apt install -y python3-pygraphviz python3-pandocfilters

clean:
	rm -f -r ../docs/graphviz-images/
	rm -f $(doc_DATA)

