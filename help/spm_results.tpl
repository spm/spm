<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
                "http://www.w3.org/TR/REC-html40/loose.dtd">
<html>
<head>
  <title>{SPM}: {CON_TITLE}</title>
  <meta name="description" content="{CON_TITLE}">
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
  <meta name="generator" content="{SPM} &copy; 1991, 1994-2013 Wellcome Trust Centre for NeuroImaging">
  <!-- <link type="text/css" rel="stylesheet" href="spm.css"> -->
  <script type="text/javascript">
    function moveto(x,y,z) {
      P1 = 125; P2 = 269; P3 = 105; P4 = 307;
      document.getElementById('cs1').style.left = P1+y+'px';
      document.getElementById('cs1').style.top  = P2+x+'px';
      document.getElementById('cs2').style.left = P1+y+'px';
      document.getElementById('cs2').style.top  = P3-z+'px';
      document.getElementById('cs3').style.left = P4+x+'px';
      document.getElementById('cs3').style.top  = P3-z+'px';
    }
  </script>
</head>
<body>

<h1 style="text-align:center;">{CON_TITLE}</h1>

<table>
  <tr>
    <td><div style="position:relative;">
      <img src="{IMG_MIP}"/>
      <!-- BEGIN cursor -->
      <img id="{CS_ID}" src="{IMG_CURSOR}" width="8" height="11" style="position:absolute;top:{CS_Y}px;left:{CS_X}px"/>
      <!-- END cursor -->
      <div style="position:absolute;top:225px;left:240px;font-size:x-large;font-weight:bold;">SPM{{STAT_STR}}</div>
    </div></td>
    <td>
      <table>
        <tr><td>
          <img src="{IMG_CON}" width="200" height="30" border="1" style="position:relative;bottom:2px;"/>
        </td></tr>
        <tr><td>
          <img src="{IMG_X}" width="200" height="300" border="1"/>
        </td></tr>
      </table>
    </td>
  </tr>
</table>

<p><strong>Statistics: <em>{RES_TITLE}</em></strong></p>
<table border="0">
  <tr>
    <td colspan="2" align="center">set-level</td><td colspan="4" align="center">cluster-level</td><td colspan="5" align="center">peak-level</td><td rowspan="2">mm mm mm</td>
  </tr>
  <tr>
    <td>p</td><td>c</td><td>p<sub>FWE-corr</sub></td><td>p<sub>FDR-corr</sub></td><td>k<sub>E</sub></td><td>p<sub>unc</sub></td><td>p<sub>FWE-corr</sub></td><td>p<sub>FDR-corr</sub></td><td>{CON_STAT}</td><td>Z<sub>E</sub></td><td>p<sub>unc</sub></td>
  </tr>
  <!-- BEGIN resrow -->
  <tr>
    <td>{RES_COL1}</td><td>{RES_COL2}</td><td>{RES_COL3}</td><td>{RES_COL4}</td><td>{RES_COL5}</td><td>{RES_COL6}</td><td>{RES_COL7}</td><td>{RES_COL8}</td><td>{RES_COL9}</td><td>{RES_COL10}</td><td>{RES_COL11}</td><td><div onmouseover="moveto({RES_XYZ});this.style.fontWeight='bold';" onmouseout="this.style.fontWeight='normal';">{RES_COL12}</div></td></tr>
  </tr>
  <!-- END resrow -->
</table>

<p><em>{RES_STR}</em></p>
<table>
  <!-- BEGIN resftr -->
  <tr><td>{RES_FTR}</td></tr>
  <!-- END resftr -->
</table>

<hr/><address>Generated on {DATE} by <strong><a href="http://www.fil.ion.ucl.ac.uk/spm/">{SPM}</a></strong> &copy; 1991, 1994-2013 Wellcome Trust Centre for NeuroImaging</address>
</body>
</html>
