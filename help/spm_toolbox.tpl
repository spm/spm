<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
                "http://www.w3.org/TR/REC-html40/loose.dtd">
<html>
<head>
  <title>{SPM} Toolboxes</title>
  <link rel="stylesheet" type="text/css" href="{SPM_CSS}" />
</head>
<body>

  <!-- BEGIN loading -->
  <h1 style="text-align:center;">{SPM} Toolboxes</h1>
  <br/><p>Loading description of toolboxes...</p>
  <p style="margin-top:5em;"><center><img src="{IMG_LOADING}"></center></p>
  <!-- END loading -->

  <!-- BEGIN listing -->
  <h1 style="text-align:center;">{SPM} Toolboxes</h1>
  <!-- BEGIN toolbox -->
  <div class="tbx">
    <h4>
      <a href="matlab:web('{TBX_URL}','-browser')">{TBX_NAME}</a>
      <span class="star"><a href="matlab:spm_toolbox('install','{TBX_ID}')"><img src="{IMG_STAR}" width="16" height="18" title="{IMG_STAR_TITLE}"/></a></span>
    </h4>
    <p><a href="mailto:{TBX_EMAIL}" title="email {TBX_AUTHOR} <{TBX_EMAIL}>">{TBX_AUTHOR}</a></p>
    <p>{TBX_SUMMARY}</p>
  </div>
  <!-- END toolbox -->

  <p>
    Found {TBX_NB} toolboxes compatible with {SPM}. 
    <span class="star"><a href="matlab:spm_toolbox;"><img src="{IMG_REFRESH}" width="16" height="18" title="Refresh"/></a></span>
  </p>
  <!-- END listing -->

</body>
</html>
