import hashlib


def hash_leaf(value):
    # Convert the value to a string and encode it as bytes
    value_str = str(value).encode("utf-8")
    # Compute the SHA-256 hash
    hash_result = hashlib.sha256(value_str).digest()
    return hash_result


def hash_children(children):
    # Create a new SHA-256 hasher
    hasher = hashlib.sha256()
    # Update the hasher with each child hash
    for child in children:
        hasher.update(child)
    # Compute the final hash
    hash_result = hasher.digest()
    return hash_result


def str_rust_u8(value):
    return "[" + ", ".join(f"{byte}" for byte in list(value)) + "]"


def calculate_tree_root(leafs, abi=False) -> bytes:
    values = [hash_leaf(str(leaf)) for leaf in leafs]
    print(f"leafs: {values}")
    print(f"first leaf: {str_rust_u8(values[0])}")
    print(f"fifth leaf: {str_rust_u8(values[5])}")

    h0 = hash_children([values[0], values[1]])
    h1 = hash_children([values[2], values[3]])
    h2 = hash_children([values[4], values[5]])
    h3 = hash_children([values[6], values[7]])
    print(
        f"1st level: {str_rust_u8(h0)}, {str_rust_u8(h1)}, {str_rust_u8(h2)}, {str_rust_u8(h3)}"
    )

    h00 = hash_children([h0, h1])
    h01 = hash_children([h2, h3])
    print(f"2nd level: {str_rust_u8(h00)}, {str_rust_u8(h01)}")

    h000 = hash_children([h00, h01])
    print(f"root: {h000}")
    u8_root = "[" + ", ".join(f"{byte}" for byte in list(h000)) + "]"
    print(f"root array: {u8_root}")

    return h000


if __name__ == "__main__":
    accounts = [x for x in range(8)]
    print(f"values: {accounts}")

    account_root = calculate_tree_root(accounts)
