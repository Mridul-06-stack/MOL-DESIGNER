import pytest
import sys

if __name__ == "__main__":
    with open("err.txt", "w") as f:
        sys.exit(pytest.main(["tests/test_scoring.py", "-k", "test_worst_case_scoring", "-v"], stdout=f))

