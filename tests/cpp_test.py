import os
import subprocess
import unittest

import pytest

SVMTK_test = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "bin", "SVMTK_test")
)


class MainTest(unittest.TestCase):
    @pytest.mark.skipif(not os.path.exists(SVMTK_test), reason="C++ tests not built")
    def test_cpp(self):
        print("\n\nTesting C++ code...")
        subprocess.check_call(SVMTK_test, shell=True)
        print()  # for prettier output


if __name__ == "__main__":
    unittest.main()
