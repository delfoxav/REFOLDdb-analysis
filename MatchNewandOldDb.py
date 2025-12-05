import pandas as pd
import sqlite3

# Load the downloaded DB as a pandas DF
new_db = pd.read_table("REFOLD_161028.txt")
print(new_db["date"].head())

# Connect to the sqlite
conn = sqlite3.connect('REFOLD.db')

old_db = pd.read_sql("SELECT * FROM main", conn)
print(old_db.head())