from pyproj import Transformer

# Criar o objeto transformador para converter de WGS 84 para Lambert-93
transformer = Transformer.from_crs(4326, 2154, always_xy=True)

# Coordenadas geográficas (longitude, latitude) de Paris
longitude = 3
latitude = 46.5

# Transformar para o sistema Lambert-93
x, y = transformer.transform(longitude, latitude)

print(f"Coordenadas projetadas em Lambert-93: x={x}, y={y}")