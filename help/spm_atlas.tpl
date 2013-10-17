<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
                "http://www.w3.org/TR/REC-html40/loose.dtd">
<html>
<head>
  <title>{SPM} Atlases</title>
  <link rel="stylesheet" type="text/css" href="{SPM_CSS}" />
</head>
<body>

  <!-- BEGIN loading -->
  <h1 style="text-align:center;">{SPM} Atlases</h1>
  <br/><p>{TXT_LOADING}</p>
  <p style="margin-top:5em;"><center><img src="{IMG_LOADING}"></center></p>
  <!-- END loading -->

  <!-- BEGIN listing -->
  <h1 style="text-align:center;">{SPM} Atlases</h1>
  <!-- BEGIN atlas -->
  <div class="tbx">
    <h3>
      <a href="matlab:web('{ATLAS_URL}','-browser');">{ATLAS_NAME}</a>
    </h3>
    <p><a href="mailto:{ATLAS_EMAIL}" title="email {ATLAS_AUTHOR} <{ATLAS_EMAIL}>">{ATLAS_AUTHOR}</a></p>
    <p>{ATLAS_SUMMARY}</p>
    <p><center><a href="matlab:spm_atlas('install','{ATLAS_ID}');">Install</a></center></p>
  </div>
  <!-- END atlas -->

  <!-- END listing -->

</body>
</html>
