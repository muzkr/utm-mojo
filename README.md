# utm-mojo

A Mojo implementation of UTM (Universal Transverse Mercator) projection for converting between geographic coordinates (latitude/longitude) and UTM grid coordinates.

**Features**: 

- **WGS84 ellipsoid support** by default with configurable ellipsoid parameters
- **Bidirectional conversion**:
  - `geo_to_grid`: Latitude/longitude -> UTM (zone, northing, easting, convergence, scale)
  - `grid_to_geo`: UTM -> Latitude/longitude
- **Simplified variants** that omit convergence and scale factor for faster computation
- **Norway/Svalbard zone exceptions**: Special UTM zone handling for high-latitude regions
- **Precision tuning**: Float32 uses 10mm tolerance, Float64 uses 1mm tolerance for inverse conversion
- **Zero dependencies**: Pure Mojo implementation

## Setup

Create a Python virtual environment and install Mojo:

```bash
# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install Mojo
pip install -r requirements.txt
```

## Testing

Test cases are provided in `test/` to verify implementation correctness and demonstrate API usage. 

Run tests:

```bash
mojo -I src test/test_utm.mojo
```

## License

MIT License - see [LICENSE](LICENSE) file

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
