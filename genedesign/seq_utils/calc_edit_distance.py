def calculate_edit_distance(s1, s2):
    """
    Compute the edit distance between two strings using a dynamic programming approach based on the Smith-Waterman algorithm for local alignment.

    Parameters:
        s1 (str): The first string to compare.
        s2 (str): The second string to compare.

    Returns:
        int: The edit distance between the two strings, defined as the minimum number of edits (insertions, deletions, or substitutions) required to transform one string into the other.
    """
    s1_len = len(s1)
    s2_len = len(s2)
    dist = [[0] * (s2_len + 1) for _ in range(s1_len + 1)]

    # Initialize distances for transformations involving empty strings
    for i in range(s1_len + 1):
        dist[i][0] = i
    for j in range(s2_len + 1):
        dist[0][j] = j

    # Compute distances
    for i in range(1, s1_len + 1):
        for j in range(1, s2_len + 1):
            if s1[i - 1] == s2[j - 1]:
                dist[i][j] = dist[i - 1][j - 1]
            else:
                dist[i][j] = 1 + min(dist[i - 1][j], dist[i][j - 1], dist[i - 1][j - 1])

    return dist[s1_len][s2_len]

def main():
    # Example usage
    pairs = [
        ("AACAAGATAT", "AACATGATAT", "Edit distance 1"),
        ("AACAAGTTAT", "ATCAAGTTCT", "Edit distance 2")
    ]
    
    for s1, s2, label in pairs:
        distance = calculate_edit_distance(s1, s2)
        print(f"{label}: {distance}")

if __name__ == "__main__":
    main()
