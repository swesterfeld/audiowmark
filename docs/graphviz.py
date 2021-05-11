#!/usr/bin/env python3
# Licensed BSD-3-Clause: https://spdx.org/licenses/BSD-3-Clause.html
"""
Pandoc filter to process code blocks with class "graphviz" into
SVG based graphviz-generated images.

Needs pygraphviz
"""
# Based on https://github.com/jgm/pandocfilters/blob/master/examples/graphviz.py

import sys, os, pygraphviz
from pandocfilters import toJSONFilter, Para, Image, RawBlock, get_filename4code, get_caption, get_extension, get_value

# https://pandoc.org/lua-filters.html

def graphviz (key, value, frmat, _):
  if key == 'CodeBlock':
    [[ident, classes, keyvals], code] = value
    if "graphviz" in classes:
      # sys.stderr.write ('Debug: ' + str (value) + '\n')
      caption, typef, keyvals = get_caption (keyvals)
      prog, keyvals = get_value (keyvals, u"prog", u"dot")
      force, keyvals = get_value (keyvals, u"force", u"false")

      filetype = get_extension (frmat, "svg", html ="svg", latex ="pdf")
      dest = get_filename4code ("graphviz", str (value), filetype)

      if force != u"false" or not os.path.isfile (dest):
        g = pygraphviz.AGraph (string = code)
        g.layout ()
        g.draw (dest, prog = prog)
        sys.stderr.write ('Create image: ' + dest + '\n')
      if dest.endswith ('.svg'):
        return RawBlock (u'html', '<object type="image/svg+xml" data="%s"></object>' % dest)
      return Para ([Image ([ident, [], keyvals], caption, [dest, typef])])

if __name__ == "__main__":
  toJSONFilter (graphviz)
