import geopandas as gpd
import zipfile
from shapely.geometry import Polygon, LineString, Point, box, MultiPolygon, MultiLineString, LinearRing
from shapely.errors import ShapelyError
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
import re


def read_kmz(filename):
    with zipfile.ZipFile(filename, 'r') as kmz:
        # Assuming only one KML file in the KMZ, if multiple, adjust accordingly
        kml_file = kmz.namelist()[0]
        with kmz.open(kml_file) as kml:
            kmldata = kml.read()

    return kmldata


def parse_kmz(kmz_data, filter_types=None):
    """
    Parses a KMZ file and extracts features based on specified geometry types.
    
    Parameters:
    kmz_data (str): The KMZ data as a string.
    filter_types (str or list of str, optional): The types of geometries to include in the output.
                                                 Can be 'Point', 'Polygon', 'LineString' (case insensitive) 
                                                 or a list of these types. If None, all types are returned.
    
    Returns:
    list of dict: A list of dictionaries where each dictionary represents a feature with 
                  geometry and properties.
    
    Example:
    >>> kmz_data = '<kml><Placemark>...</Placemark></kml>'
    >>> parse_kmz(kmz_data, filter_types='Polygon')
    [{'geometry': <shapely.geometry.polygon.Polygon object at 0x...>, 'name': '...', 'description': '...'}]
    >>> parse_kmz(kmz_data, filter_types=['Polygon', 'Point'])
    [{'geometry': <shapely.geometry.polygon.Polygon object at 0x...>, 'name': '...', 'description': '...'}, 
     {'geometry': <shapely.geometry.point.Point object at 0x...>, 'name': '...', 'description': '...'}]
    """
    features = []
    tree = ET.ElementTree(ET.fromstring(kmz_data))
    root = tree.getroot()

    if filter_types:
        if isinstance(filter_types, str):
            filter_types = [filter_types]
        filter_types = [ft.lower() for ft in filter_types]
    else:
        filter_types = ['point', 'polygon', 'linestring']

    for placemark in root.findall(".//{http://www.opengis.net/kml/2.2}Placemark"):
        geom = None
        properties = {}
        for element in placemark.iter():
            if element.tag.endswith('Point') and 'point' in filter_types:
                coordinates_text = element.find(
                    './/{http://www.opengis.net/kml/2.2}coordinates').text.strip()
                lon, lat, _ = map(float, coordinates_text.split(','))
                geom = Point(lon, lat)
            elif element.tag.endswith('Polygon') and 'polygon' in filter_types:
                coordinates_text = element.find(
                    './/{http://www.opengis.net/kml/2.2}coordinates').text.strip()
                coordinates = [tuple(map(float, coord.split(',')))
                               for coord in coordinates_text.split()]
                geom = Polygon(coordinates)
            elif element.tag.endswith('LineString') and 'linestring' in filter_types:
                coordinates_text = element.find(
                    './/{http://www.opengis.net/kml/2.2}coordinates').text.strip()
                coordinates = [tuple(map(float, coord.split(',')))
                               for coord in coordinates_text.split()]
                geom = LineString(coordinates)
            elif element.tag.endswith('name'):
                properties['name'] = element.text.strip()
            elif element.tag.endswith('description'):
                properties['description'] = element.text.strip()

        if geom is not None:
            properties['geometry'] = geom
            features.append(properties)
            
    return gpd.GeoDataFrame(features)



def extract_links_from_html(html_content):
    """
    Extract links from HTML content and return a dictionary.

    Parameters:
    html_content (str): HTML content as a string.

    Returns:
    dict: A dictionary where keys are strings in <th> tags and values are lists of hrefs in the corresponding <td> tags.
    """
    soup = BeautifulSoup(html_content, 'html.parser')
    result = {}

    # Find all rows in the table
    for tr in soup.find_all('tr'):
        th = tr.find('th')
        if th:
            key = th.text.strip()
            result[key] = []
            for td in tr.find_all('td'):
                a_tags = td.find_all('a', href=True)
                for a in a_tags:
                    result[key].append(a['href'])
    # Remove entries with empty lists
    result = {k: v for k, v in result.items() if v}

    return result


def filter_geodataframe(gdf, bbox=None, shapefile_path=None, remove_date_line_crossing=True):
    """
    Filter a GeoDataFrame based on a bounding box or a shapefile and optionally remove geometries crossing the International Date Line.

    Parameters:
    gdf (geopandas.GeoDataFrame): The GeoDataFrame to filter.
    bbox (list, optional): A bounding box specified as [min_lat, min_lon, max_lat, max_lon]. Either bbox or shapefile_path must be provided.
    shapefile_path (str, optional): Path to a shapefile for filtering. Either bbox or shapefile_path must be provided.
    remove_date_line_crossing (bool, optional): If True, remove geometries crossing the International Date Line. Default is True.

    Returns:
    geopandas.GeoDataFrame: The filtered GeoDataFrame.
    """
    if (bbox and shapefile_path) or (not bbox and not shapefile_path):
        raise ValueError("Either bbox or shapefile_path must be provided, but not both")

    try:
        if bbox:
            # Ensure bbox is in the format [min_lat, min_lon, max_lat, max_lon]
            if len(bbox) != 4:
                raise ValueError("Bounding box must be a list of 4 elements: [min_lat, min_lon, max_lat, max_lon]")

            min_lon, min_lat, max_lon, max_lat = bbox
            # Create a bounding box from the given coordinates
            bounding_box = box(min_lon, min_lat, max_lon, max_lat)

            # Create a spatial index for the GeoDataFrame
            sindex = gdf.sindex

            # Find all geometries in gdf that intersect the bounding box
            possible_matches_index = list(sindex.intersection(bounding_box.bounds))
            possible_matches = gdf.iloc[possible_matches_index]

            # Perform a precise intersection
            filtered_gdf = possible_matches[possible_matches.intersects(bounding_box)]

        elif shapefile_path:
            try:
                # Read the shapefile
                area_gdf = gpd.read_file(shapefile_path)
            except Exception as e:
                raise ValueError(f"Error reading shapefile: {e}")

            # Ensure the shapefile contains only one geometry (if multiple, combine them)
            if len(area_gdf) > 1:
                area_gdf = area_gdf.unary_union
            else:
                area_gdf = area_gdf.geometry.iloc[0]

            # Create a spatial index for the GeoDataFrame
            sindex = gdf.sindex

            # Find all geometries in gdf that intersect the shapefile geometry
            possible_matches_index = list(sindex.intersection(area_gdf.bounds))
            possible_matches = gdf.iloc[possible_matches_index]

            # Perform a precise intersection
            filtered_gdf = possible_matches[possible_matches.intersects(area_gdf)]

        # Optionally remove geometries crossing the International Date Line
        if remove_date_line_crossing:
            filtered_gdf = filtered_gdf[~filtered_gdf['geometry'].apply(has_horizontal_portion_crossing_date_line)]

        return filtered_gdf

    except ValueError as ve:
        raise ve
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return gpd.GeoDataFrame()  # Return an empty GeoDataFrame on error


def has_horizontal_portion_crossing_date_line(geometry):
    """
    Determine if a geometric shape has a horizontal portion crossing the International Date Line (longitude ±180 degrees).

    Parameters:
    geometry (shapely.geometry.base.BaseGeometry): A Shapely geometry object which can be a LineString, Polygon, 
                                                   MultiLineString, or MultiPolygon.

    Returns:
    bool: True if the geometry crosses the International Date Line, False otherwise.
    """
    def extract_lon_lat(coord):
        """Extract longitude and latitude from coordinate tuple, handling 2D and 3D coordinates."""
        if len(coord) == 2:
            return coord
        elif len(coord) == 3:
            return coord[0], coord[1]
        else:
            raise ValueError("Invalid coordinate dimension")

    try:
        if isinstance(geometry, LineString) or isinstance(geometry, LinearRing):
            coords = list(geometry.coords)
            lon1, lat1 = extract_lon_lat(coords[0])
            lon2, lat2 = extract_lon_lat(coords[-1])
            # Check if the LineString crosses the International Date Line
            if abs(lon1 - lon2) > 180:
                return True
            else:
                return False
        elif isinstance(geometry, Polygon):
            # Check the exterior ring of the Polygon
            if has_horizontal_portion_crossing_date_line(geometry.exterior):
                return True
            # Check all interior rings (holes) of the Polygon
            for interior in geometry.interiors:
                if has_horizontal_portion_crossing_date_line(interior):
                    return True
            return False
        elif isinstance(geometry, MultiLineString):
            # Check each LineString in the MultiLineString
            for part in geometry.geoms:
                if has_horizontal_portion_crossing_date_line(part):
                    return True
            return False
        elif isinstance(geometry, MultiPolygon):
            # Check each Polygon in the MultiPolygon
            for part in geometry.geoms:
                if has_horizontal_portion_crossing_date_line(part):
                    return True
            return False
        else:
            # If the geometry is not one of the types checked above, return False
            return False
    except (AttributeError, ShapelyError, ValueError) as e:
        # Handle errors related to attribute access, Shapely operations, and invalid coordinates
        print(f"Error processing geometry: {e}")
        return False




def extract_pass_scene_numbers(html, name):
    pass_match = re.search(r'Pass:</b> (\d+)', html)
    scene_match = re.search(r'Scene:</b> (\d+[A-Z])', html)
    
    if pass_match:
        pass_number = pass_match.group(1)
    else:
        # Try to extract pass number from the name if not found in the description
        pass_number_match = re.search(r'Pass (\d+)', name)
        pass_number = pass_number_match.group(1) if pass_number_match else None
        
    scene_number = scene_match.group(1) if scene_match else None
    return pass_number, scene_number

def process_dataframe(gdf):
    # Initialize lists to store extracted pass and scene numbers
    pass_numbers = []
    scene_numbers = []

    # Loop through each description to extract pass and scene numbers
    for i, row in gdf.iterrows():
        pass_number, scene_number = extract_pass_scene_numbers(row['description'], row['name'])
        pass_numbers.append(pass_number)
        scene_numbers.append(scene_number)

    # Add the results as new columns to the dataframe
    gdf['pass_number'] = pass_numbers
    gdf['scene_number'] = scene_numbers

    # Drop the original 'description' column if no longer needed
    gdf.drop(columns=['description'], inplace=True)

    return gdf