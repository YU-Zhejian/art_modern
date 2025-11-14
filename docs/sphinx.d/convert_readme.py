import re
import os
import hashlib
import time

import requests

IMG_REGEX = re.compile(r"\[!\[(.+)\]\((.+)\)\]\((.+)\)")
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)

fetch_lru = {}

# Deal with Readme.md
with open(os.path.join(THIS_DIR, "src", "Readme.md"), encoding="UTF-8") as f:
    with open(os.path.join(THIS_DIR, "src", "README_NEW.md"), "w", encoding="UTF-8") as w:
        for l in f:
            # Match for lines like [![Anaconda-Server Badge](https://anaconda.org/bioconda/art_modern/badges/downloads.svg)](https://anaconda.org/bioconda/art_modern)
            m = IMG_REGEX.match(l)
            if m:
                alt_text = m.group(1)
                img_url = m.group(2)
                link_url = m.group(3)
                if img_url in fetch_lru:
                    w.write(f"[![{alt_text}]({fetch_lru[img_url]})]({link_url})\n")
                    continue
                # Try to fetch the image to see if it's valid
                try:
                    response = requests.get(img_url, timeout=60)
                    print( f"Fetched image URL: {img_url}" )
                    if response.status_code == 200:
                        # Get the SHA256 hash of the image content
                        img_data = response.content
                        # Tell whether the data is SVG
                        is_svg = img_data.lstrip().startswith(b"<svg")
                        if is_svg:
                            digest = hashlib.sha256(img_data).hexdigest()
                            if b"Unable to select next GitHub token from pool" in img_data:
                                print( f"Failed to fetch image URL: {img_url} with status code {response.status_code}" )
                                w.write(f"[{alt_text}]({link_url})\n")
                                continue
                            new_file_name = f"{digest}.svg"
                            fetch_lru[img_url] = new_file_name
                            new_path = os.path.join(THIS_DIR, "src", new_file_name)
                            with open(new_path, "wb") as img_file:
                                img_file.write(img_data)
                            # Write the line with the original image URL
                            w.write(f"[![{alt_text}]({new_file_name})]({link_url})\n")
                        else:
                            # For non-SVG images, just write the original line
                            w.write(f"[{alt_text}]({link_url})\n")
                    else:
                        print( f"Failed to fetch image URL: {img_url} with status code {response.status_code}" )
                        w.write(f"[{alt_text}]({link_url})\n")
                except requests.RequestException as e:
                    print( f"Error fetching image URL: {img_url}. Error: {e}" )
                    raise e
                time.sleep(10)
            else:
                w.write(l)
