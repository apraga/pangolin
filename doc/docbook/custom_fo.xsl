<?xml version='1.0'?> 
<xsl:stylesheet  
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  version="1.0"> 

  <!-- Split into different HTML -->
  <xsl:import href="/usr/share/sgml/docbook/xsl-stylesheets/fo/docbook.xsl"/> 

  <!-- Correct command style -->
  <xsl:template match="command">
    <xsl:call-template name="inline.monoseq"/>
  </xsl:template>

  <!-- TOC after abstract -->
  <xsl:param name="generate.toc"/>
  <xsl:param name="process.empty.source.toc" select="1"/>
</xsl:stylesheet>  
