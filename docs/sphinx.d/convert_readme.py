import io
import re
import os
import hashlib
import time
from typing import Mapping

import requests

IMG_REGEX = re.compile(r"\[!\[(.+)\]\((.+)\)\]\((.+)\)")
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)


def get_file(
    w: io.TextIOBase,
    img_url: str,
    alt_text: str,
    link_url: str,
    fetch_lru: Mapping[str, str],
):
    try:
        response = requests.get(img_url, timeout=60)
        print(f"Fetched image URL: {img_url}")
        if response.status_code == 200:
            # Get the SHA256 hash of the image content
            img_data = response.content
            # Tell whether the data is SVG
            is_svg = img_data.lstrip().startswith(b"<svg")
            if is_svg:
                digest = hashlib.sha256(img_data).hexdigest()
                if b"Unable to select next GitHub token from pool" in img_data:
                    print(
                        f"Failed to fetch image URL: {img_url} with status code {response.status_code}: Upstream service failed."
                    )
                    w.write(f"[{alt_text}]({link_url})\n")
                    return
                new_file_name = f"{digest}.svg"
                fetch_lru[img_url] = new_file_name
                new_path = os.path.join(THIS_DIR, "src", new_file_name)
                with open(new_path, "wb") as img_file:
                    img_file.write(img_data)
                # Write the line with the original image URL
                w.write(f"[![{alt_text}]({new_file_name})]({link_url})\n")
                # Wait on each successful fetch to avoid rate limiting
                time.sleep(10)
            else:
                # For non-SVG images, just write the original line
                w.write(f"[![{alt_text}]({img_url})]({link_url})\n")
        else:
            print(f"Failed to fetch image URL: {img_url} with status code {response.status_code}")
            w.write(f"[{alt_text}]({link_url})\n")
    except requests.RequestException as e:
        print(f"Error fetching image URL: {img_url}. Error: {e}")
        w.write(f"[{alt_text}]({link_url})\n")


if __name__ == "__main__":
    fetch_lru_ = {}
    release = os.environ.get("PACKAGE_VERSION")
    if not release:
        release = "0.0.0-dev"

    CACHE_RELEASE_FILE = os.path.join(THIS_DIR, f"fetch_lru_{release}.txt")
    if os.path.exists(CACHE_RELEASE_FILE):
        with open(CACHE_RELEASE_FILE, "r", encoding="UTF-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    img_url_, file_name_ = parts
                    fetch_lru_[img_url_] = file_name_

    # Deal with Readme.md
    with open(os.path.join(THIS_DIR, "src", "Readme.md"), encoding="UTF-8") as f:
        with open(os.path.join(THIS_DIR, "src", "README_NEW.md"), "w", encoding="UTF-8") as w:
            for l in f:
                # Match for lines like [![Anaconda-Server Badge](https://anaconda.org/bioconda/art_modern/badges/downloads.svg)](https://anaconda.org/bioconda/art_modern)
                m = IMG_REGEX.match(l)
                if m:
                    alt_text_ = m.group(1)
                    img_url_ = m.group(2)
                    link_url_ = m.group(3)
                    if img_url_ in fetch_lru_:
                        file_name_ = fetch_lru_[img_url_]
                        if os.path.exists(os.path.join(THIS_DIR, "src", file_name_)):
                            # Write the line with the cached image file
                            w.write(f"[![{alt_text_}]({fetch_lru_[img_url_]})]({link_url_})\n")
                            continue
                        else:
                            print(f"LRU Cache miss for image URL: {img_url_}, refetching.")
                    # Try to fetch the image to see if it's valid
                    get_file(w, img_url_, alt_text_, link_url_, fetch_lru_)
                else:
                    w.write(l)
    # Write back the updated fetch_lru_
    with open(CACHE_RELEASE_FILE, "w", encoding="UTF-8") as f:
        for img_url_, file_name_ in fetch_lru_.items():
            f.write(f"{img_url_}\t{file_name_}\n")
