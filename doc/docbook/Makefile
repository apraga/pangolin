STYLE_HTML=custom_html.xsl  
STYLE_FO=custom_fo.xsl  

# Select the location of the stylesheet file here
LOC_HTML=/usr/share/sgml/docbook/xsl-stylesheets/html/chunk.xsl
#LOC_HTML=/usr/share/xml/docbook/xsl-stylesheets-1.78.1/html/chunk.xsl
#LOC_HTML=/usr/share/xml/docbook/stylesheet/docbook-xsl/html/chunk.xsl

# Change stylesheet for fop
TMP=$(subst html,fo,${LOC_HTML})
# Chank form chunk to docbook
LOC_FO=$(subst chunk,docbook,${TMP})

#Escape backslash for subsitution
LOC_HTML2=$(subst /,\/,${LOC_HTML})
LOC_FO2=$(subst /,\/,${LOC_FO})

# Use FreeBSD CSS
OPTS=--stringparam html.stylesheet docbook_freebsd.css
DEST=user_manual

all: html
	
html: replace_html ${DEST}/index.html

pdf: replace_fo pangolin_manual.pdf

pangolin_manual.pdf:
	xsltproc ${LOC_FO} pangolin-doc.xml > tmp.fo
	fop tmp.fo $@

${DEST}/index.html: pangolin-doc.xml
	xsltproc ${OPTS} -o ${DEST}/ ${STYLE_HTML} $?

# Replace stylesheet directory in  our custom stylesheet files
replace_html:
	perl -pi -e 's/(?<=xsl:import href=").*(?="\/>)/${LOC_HTML2}/' ${STYLE_HTML}

replace_fo:
	perl -pi -e 's/(?<=xsl:import href=").*(?="\/>)/${LOC_FO2}/' ${STYLE_FO}

.PHONY: ${DEST}/index.html
.PHONY: pangolin_manual.pdf

