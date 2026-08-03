dist_man_MANS = audiowmark.1 videowmark.1

rebuild-manpages:
	for p in $(dist_man_MANS); do \
          pandoc $$p.md -s -M date="$$(date +%Y-%m-%d)" \
                           -M section=1 \
                           -t man > $$p; \
        done
