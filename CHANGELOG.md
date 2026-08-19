# SunCET Instrument Simulator Code Base Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Added
* Nothing yet

### Changed
* Model 16-bit firmware overflow as modulo wraparound instead of saturation,
  while retaining 32-bit particle-filter accumulation before the configured
  right shift back into the output buffer.

### Deprecated
* Nothing yet

### Removed
* Nothing yet

### Fixed
* Nothing yet
