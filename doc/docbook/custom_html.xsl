<?xml version='1.0'?> 
<xsl:stylesheet  
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  version="1.0"> 

  <!-- Split into different HTML -->
  <xsl:import href="/usr/share/sgml/docbook/xsl-stylesheets/html/chunk.xsl"/> 

  <!-- Correct command style -->
  <xsl:template match="command">
    <xsl:call-template name="inline.monoseq"/>
  </xsl:template>

  <!-- TOC after abstract -->
  <xsl:param name="generate.toc"/>
  <xsl:param name="process.empty.source.toc" select="1"/>

  <!-- No break after table -->
  <xsl:param name="formal.object.break.after">0</xsl:param>
</xsl:stylesheet>  
