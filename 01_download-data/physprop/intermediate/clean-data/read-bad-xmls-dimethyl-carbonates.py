import pandas as pd
from xml.etree import ElementTree
import re
import tqdm

# First extract all variable definitions per PureOrMixtureData block

def read_xml(file):
    with open(file, "r") as f:
        xml = f.read()

    linked_entries = []
    
    root = ElementTree.fromstring(xml)
    
    namespace_string = re.search(r"{.*\}", root.tag).group(0)[1:-1]
    ns = {"ThermoML": namespace_string}

    compound_map = {}
    for compound in root.findall('ThermoML:Compound', ns):
        org_num = compound.find('ThermoML:RegNum/ThermoML:nOrgNum', ns)
        names = compound.findall('ThermoML:sCommonName', ns)
        if org_num is not None and names:
            compound_map[int(org_num.text)] = names[0].text
        
    for pod in root.findall('ThermoML:PureOrMixtureData', ns):
        # Build variable map: var_number -> (compound name, phase, type)
        variable_map = {}
        for var in pod.findall('ThermoML:Variable', ns):
            var_number = int(var.find('ThermoML:nVarNumber', ns).text)
            try:
                org_num = int(var.find('ThermoML:VariableID/ThermoML:RegNum/ThermoML:nOrgNum', ns).text)
            except:
                continue
            compound_name = compound_map.get(org_num, f"Unknown {org_num}")
            var_type = var.find('ThermoML:VariableID/ThermoML:VariableType/ThermoML:eComponentComposition', ns)
            var_type_text = var_type.text if var_type is not None else "Unknown type"
            variable_map[var_number] = {
                "compound": compound_name,
                "type": var_type_text
            }
    
        # Get substances (for reference only)
        substances = []
        for comp in pod.findall('ThermoML:Component', ns):
            org_num = int(comp.find('ThermoML:RegNum/ThermoML:nOrgNum', ns).text)
            substances.append(compound_map.get(org_num, f"Unknown {org_num}"))
    
        # Get property type
        prop = pod.find('ThermoML:Property', ns)
        prop_type_elem = prop.find('.//ThermoML:ePropName', ns)
        prop_type = prop_type_elem.text if prop_type_elem is not None else 'Unknown'
    
        # Link variable value with associated variable metadata
        for nv in pod.findall('ThermoML:NumValues', ns):
            val_elem = nv.find('ThermoML:PropertyValue/ThermoML:nPropValue', ns)
            if val_elem is None:
                continue
            value = float(val_elem.text)
    
            for var_val in nv.findall('ThermoML:VariableValue', ns):
                var_number = int(var_val.find('ThermoML:nVarNumber', ns).text)
                mole_fraction = float(var_val.find('ThermoML:nVarValue', ns).text)
                var_info = variable_map.get(var_number, {})
                linked_entries.append({
                    "Property type": prop_type,
                    "Value": value,
                    "Substances": substances,
                    "Mole fraction compound": var_info.get("compound", "Unknown"),
                    "Mole fraction type": var_info.get("type", "Unknown"),
                    "Mole fraction value": mole_fraction
                })
    
    df_linked = pd.DataFrame(linked_entries)
    return df_linked

if __name__ == "__main__":
    BAD_DOIS = [
        "10.1021/je050052+",
    ]
    all_dfs = []
    for doi in BAD_DOIS:
        df = read_xml(f"{doi}.xml")
        df["doi"] = doi
        all_dfs.append(df)


    df = pd.concat(all_dfs)
    # in this paper, mole fractions should be specified with dimethyl carbonate but some have methanol
    bad = df[df["Mole fraction compound"] == "methanol"]
    bad.to_csv("bad-dimethyl-carbonate-data.csv")
