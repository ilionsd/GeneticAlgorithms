#!/usr/bin/env bash

ROOT=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JSONCPP_VERSION="1.9.5"
DOWNLOAD_DIR="$ROOT/download"
SOURCES_DIR="$DOWNLOAD_DIR/jsoncpp/$JSONCPP_VERSION"

REMOTE="https://github.com/open-source-parsers/jsoncpp/archive/refs/tags/$JSONCPP_VERSION.zip"
LOCAL="$DOWNLOAD_DIR/jsoncpp-$JSONCPP_VERSION.zip"

if [ -f "$LOCAL" ]; then
    echo "Skipping download"
else
    mkdir -p "$DOWNLOAD_DIR"
    curl -o "$LOCAL" -L "$REMOTE"
fi

if [ -d "$SOURCES_DIR" ]; then
    rm -rf "$SOURCES_DIR"
fi
mkdir -p "$SOURCES_DIR"
unzip -q "$LOCAL" "jsoncpp-$JSONCPP_VERSION/*" -d "$SOURCES_DIR"
pushd "$SOURCES_DIR/jsoncpp-$JSONCPP_VERSION"
mv * ../
popd
rm -rf "$SOURCES_DIR/jsoncpp-$JSONCPP_VERSION"

INSTALL_PREFIX="$ROOT/jsoncpp"
if [ -d "$INSTALL_PREFIX" ]; then
    rm -rf "$INSTALL_PREFIX"
fi

mkdir -p "$INSTALL_PREFIX"

BUILD_DIR="$SOURCES_DIR/build"
if [ -d "$BUILD_DIR" ]; then
    rm -rf "$BUILD_DIR"
fi
mkdir "$BUILD_DIR"

pushd "$BUILD_DIR"
cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" ../
make install
popd
