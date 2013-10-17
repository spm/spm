<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
                "http://www.w3.org/TR/REC-html40/loose.dtd">
<html>
<head>
  <title>{SPM}: {CON_TITLE}</title>
  <meta name="description" content="{CON_TITLE}">
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
  <meta name="generator" content="{SPM} &copy; 1991, 1994-{YEAR} Wellcome Trust Centre for NeuroImaging">
  <!-- <link type="text/css" rel="stylesheet" href="spm.css"> -->
  <style type="text/css">
    .crisp {
      image-rendering: -moz-crisp-edges;         /* Firefox */
      image-rendering:   -o-crisp-edges;         /* Opera */
      image-rendering: -webkit-optimize-contrast;/* Webkit (non-standard naming) */
      image-rendering: crisp-edges;
      -ms-interpolation-mode: nearest-neighbor;  /* IE (non-standard property) */
    }
  </style>
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

<!--<h1 style="text-align:center;">{CON_TITLE}</h1>-->

<table>
  <tr>
    <td colspan="2" align="center"><h1 style="text-align:center;">{CON_TITLE}</h1></td>
  </tr>
  <tr>
    <td><div style="position:relative;">
      <img src="{IMG_MIP}" class="crisp"/>
      <!-- BEGIN cursor -->
      <img id="{CS_ID}" src="{IMG_CURSOR}" class="crisp" width="8" height="11" style="position:absolute;top:{CS_Y}px;left:{CS_X}px"/>
      <!-- END cursor -->
      <div style="position:absolute;top:225px;left:240px;font-size:x-large;font-weight:bold;">SPM{{STAT_STR}}</div>
    </div></td>
    <td>
      <table>
        <tr><td>
          <img src="{IMG_CON}" class="crisp" width="200" height="30" border="1" style="position:relative;bottom:2px;"/>
        </td></tr>
        <tr><td>
          <img src="{IMG_X}" class="crisp" width="200" height="300" border="1"/>
        </td></tr>
      </table>
    </td>
  </tr>
</table>

<!--<p><strong>Statistics: <em>{RES_TITLE}</em></strong></p>-->
<table border="0">
  <tr>
    <td colspan="12" style="border-bottom: 5px solid #F00;"><strong>Statistics: <em>{RES_TITLE}</em></strong></td>
  </tr>
  <tr>
    <td colspan="2" align="center" style="border-bottom: 2px solid #F00;">set-level</td>
    <td colspan="4" align="center" style="border-bottom: 2px solid #F00;">cluster-level</td>
    <td colspan="5" align="center" style="border-bottom: 2px solid #F00;">peak-level</td>
    <td rowspan="2" align="center" style="padding-left:5px;padding-right:5px;">mm mm mm</td>
  </tr>
  <tr>
    <td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em></td><td align="center" style="padding-left:5px;padding-right:5px;">c</td>
    <td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em><sub>FWE-corr</sub></td><td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em><sub>FDR-corr</sub></td><td align="center" style="padding-left:5px;padding-right:5px;">k<sub>E</sub></td><td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em><sub>unc</sub></td>
    <td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em><sub>FWE-corr</sub></td><td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em><sub>FDR-corr</sub></td><td align="center" style="padding-left:5px;padding-right:5px;">{CON_STAT}</td><td align="center" style="padding-left:5px;padding-right:5px;">Z<sub>E</sub></td><td align="center" style="padding-left:5px;padding-right:5px;"><em>p</em><sub>unc</sub></td>
  </tr>
  <tr>
    <td colspan="12" style="border-bottom: 2px solid #F00;"><span style="font-size:1px;">&nbsp;</span></td>
  </tr>
  <!-- BEGIN resrow -->
  <tr>
    <td align="center">{RES_COL1}</td><td align="center">{RES_COL2}</td><td align="center">{RES_COL3}</td><td align="center">{RES_COL4}</td><td align="center">{RES_COL5}</td><td align="center">{RES_COL6}</td><td align="center">{RES_COL7}</td><td align="center">{RES_COL8}</td><td align="center">{RES_COL9}</td><td align="center">{RES_COL10}</td><td align="center">{RES_COL11}</td><td align="center"><div onmouseover="moveto({RES_XYZ});this.style.fontWeight='bold';" onmouseout="this.style.fontWeight='normal';">{RES_COL12}</div></td></tr>
  </tr>
  <!-- END resrow -->
</table>

<p><em>{RES_STR}</em></p>
<table>
  <!-- BEGIN resftr -->
  <tr><td>{RES_FTR}</td></tr>
  <!-- END resftr -->
</table>

<hr/><address>Generated on {DATE} by <strong><a href="http://www.fil.ion.ucl.ac.uk/spm/">{SPM}</a></strong> &copy; 1991, 1994-{YEAR} Wellcome Trust Centre for NeuroImaging</address>
</body>
</html>
