import urllib.request
import json
import sys

def fetch_and_convert_to_tremer(haplogroup_url):
    print("Fetching data from FTDNA...")
    
    try:
        # 1. Download the JSON data from the FTDNA API
        req = urllib.request.Request(haplogroup_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            
        print("\n/SNPTREE")
        
        # 2. Extract the main haplogroup node
        name = data.get('name', '')
        
        # FTDNA stores age data differently depending on the node, 
        # usually under 'age', 'tmrca', or 'formation'
        age = 0
        if 'tmrca' in data and data['tmrca']:
            # Sometimes stored as a dictionary with an 'estimate'
            age = data['tmrca'].get('estimate', 0) if isinstance(data['tmrca'], dict) else data.get('tmrca', 0)
        
        # Extract the parent
        parent = data.get('parent', {}).get('name', '') if isinstance(data.get('parent'), dict) else data.get('parent', '')

        # Clean the "R-" prefix if it exists (TREMER handles base names better)
        if name.startswith("R-"): name = name[2:]
        if parent.startswith("R-"): parent = parent[2:]

        # 3. Print the TREMER formatted line
        if name:
            print(f"{name} {age} {parent}")
            
        # 4. Optional: If the JSON includes a 'path' array (ancestors), print them too!
        if 'path' in data:
            for ancestor in data['path']:
                anc_name = ancestor.get('name', '')
                if anc_name.startswith("R-"): anc_name = anc_name[2:]
                
                anc_age = 0
                if 'tmrca' in ancestor and ancestor['tmrca']:
                    anc_age = ancestor['tmrca'].get('estimate', 0) if isinstance(ancestor['tmrca'], dict) else ancestor.get('tmrca', 0)
                
                anc_parent = ancestor.get('parent', '')
                if anc_parent.startswith("R-"): anc_parent = anc_parent[2:]
                
                if anc_name:
                    print(f"{anc_name} {anc_age} {anc_parent}")
                    
    except Exception as e:
        print(f"Error fetching or parsing data: {e}")
        sys.exit(1)

# Run the script with your URL
if __name__ == "__main__":
    url = "https://discover.familytreedna.com/resources/y-dna/R-BY17509.json"
    fetch_and_convert_to_tremer(url)
