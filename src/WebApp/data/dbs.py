# Insert your databases here with their relevant backends and passwords and stuff
{
  'ZINC':  {
    'ENGINE': 'django.db.backends.postgresql_psycopg2',
    'NAME': 'ZINCDB',
    'USER': 'postgres',
    'PASSWORD': 'RDKitWanker2007!',
    'HOST': 'rommel',
    'PORT': '5432',
    'TYPE': 'RDKit',
            },
  'MMPDB': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': r"W:\Informatics\Pharmacophore\anthony\DPhil\CODE\CHOC\src\WebApp\mmpoakley.db",#mmprealbig.db",
        # The following settings are not used with sqlite3:
        'USER': '',
        'PASSWORD': '',
        'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
        'PORT': '',                     # Set to empty string for default.
        'TYPE': 'MMP',
    }
            }
