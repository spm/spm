INSTALLATION INSTRUCTIONS FOR MACOS
===================================

1. Unzip the downloaded file (if not already unzipped).
   You should see this README file and a SPM folder.

2. Open the "Terminal" app (Command+Space, type "Terminal").

3. Paste the following command, but DO NOT press Enter yet:

   sudo xattr -cr 

   (Make sure there is a space at the end)

   [Why? This command removes the "Quarantine" flag that macOS places on downloaded files.
    Without this, macOS will block the app from running or accessing its own internal files.]

4. Drag the SPM folder into the Terminal window.
   This will automatically fill in the correct path.

5. Press Enter.
6. Type your Mac password when prompted and press Enter.

------------------------------------------------------------------

HOW TO RUN / INSTALL:

1. Install MATLAB Runtime (Required First):
   - Open the SPM folder.
   - Double-click "Runtime_Installer.app".
   - Follow the on-screen instructions.
   - Note: Please install Java (JRE), if prompted.

2. Run SPM:
   - Open the SPM folder (if not already listed).
   - Double-click "SPM.app".
   - (Optional: You can drag "SPM.app" to your Applications folder).
