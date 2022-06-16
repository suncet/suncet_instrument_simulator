import sys
from pathlib import Path
p=Path(__file__).parents[1]
sys.path.append(p.as_posix())
import example

def test_example():
    assert example.example_func(4) == 8
    assert example.example_func(3) == 6

