#!/usr/bin/env bash
# Package the compiled macOS standalone as SPM.app with the Runtime installer,
# as shipped in the release zip. Used by release.yml and install_test_standalone.yml.
#
# Usage, from the SPM checkout after spm_make_standalone:
#   package_macos_app.sh <build_tag> <zip_path>
# Reads ../standalone and ../runtime_installer.

set -ex

BUILD_TAG="$1"
ZIP_PATH="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"

# Debug: Show input artifacts
echo "Listing structure of ../standalone:"
ls -R ../standalone
# 1. Prepare Release Folder with Versioned Subfolder
APP_FOLDER_NAME="SPM-${BUILD_TAG}"
mkdir -p "release_package/SPM_MacOS_Release/$APP_FOLDER_NAME"
ROOT_DIR="release_package/SPM_MacOS_Release"
APP_DIR="$ROOT_DIR/$APP_FOLDER_NAME"

# 2. Locate and Move SPM App
SPM_APP=$(find ../standalone -maxdepth 1 -name "*.app" | head -n 1)

if [ -n "$SPM_APP" ]; then
   echo "Found SPM App: $SPM_APP"
   mv "$SPM_APP" "$APP_DIR/SPM.app"
else
   echo "Error: Could not find .app bundle in ../standalone"
   exit 1
fi

# 3. Locate and Copy MCR Installer
INSTALLER_SOURCE=$(find ../runtime_installer -name "Runtime*" | head -n 1)
if [ -n "$INSTALLER_SOURCE" ]; then
   echo "Found installer: $INSTALLER_SOURCE"
   cp -R -v "$INSTALLER_SOURCE" "$APP_DIR/Runtime_Installer.app"
else
   echo "Warning: No MCR Installer found."
fi

# 4. Icon Replacement (macOS only)
# Always download from Organization
echo "Downloading logo from 'spm' organization..."
curl -L -o "spm_org_logo.png" "https://github.com/spm.png"

if [ -s "spm_org_logo.png" ]; then
   echo "Downloaded organization logo successfully."
   ICON_SRC="spm_org_logo.png"
else
   echo "Warning: Failed to download logo."
   ICON_SRC=""
fi

if [ -f "$ICON_SRC" ]; then
   echo "Generating ICNS from $ICON_SRC..."
   mkdir spm.iconset
   # Create various sizes (sips is standard on macOS)
   sips -z 16 16     "$ICON_SRC" --out spm.iconset/icon_16x16.png
   sips -z 32 32     "$ICON_SRC" --out spm.iconset/icon_16x16@2x.png
   sips -z 32 32     "$ICON_SRC" --out spm.iconset/icon_32x32.png
   sips -z 64 64     "$ICON_SRC" --out spm.iconset/icon_32x32@2x.png
   sips -z 128 128   "$ICON_SRC" --out spm.iconset/icon_128x128.png
   sips -z 256 256   "$ICON_SRC" --out spm.iconset/icon_128x128@2x.png
   sips -z 512 512   "$ICON_SRC" --out spm.iconset/icon_512x512.png
   # Create icns
   iconutil -c icns spm.iconset -o AppIcon.icns

   # Locate Info.plist to find current icon name
   PLIST="$APP_DIR/SPM.app/Contents/Info.plist"
   # Defaults usually 'membrane.icns' or similar. We will just overwrite and update plist.
   cp AppIcon.icns "$APP_DIR/SPM.app/Contents/Resources/spm_icon.icns"

   # Update Info.plist to point to new icon
   # Using plutil to replace CFBundleIconFile
   plutil -replace CFBundleIconFile -string "spm_icon.icns" "$PLIST"
   echo "Icon updated successfully."
else
   echo "Warning: No icon source found."
fi

# 5. Copy Instructions
if [ -f "config/README_macOS_standalone.txt" ]; then
   cp "config/README_macOS_standalone.txt" "$ROOT_DIR/README_macOS_standalone.txt"
else
   echo "Warning: config/README_macOS_standalone.txt not found."
   echo "Please run: sudo xattr -cr ." > "$ROOT_DIR/README_macOS_standalone.txt"
fi

# 6. Ensure Binaries are Executable (just in case)
chmod +x "$APP_DIR/SPM.app/Contents/MacOS/"* || true

# 7. Zip Release
echo "Zipping release..."
cd "release_package"
zip -r -y "$ZIP_PATH" ./*
