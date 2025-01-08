import numpy as np
import numpy.typing as npt

def get(data: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    print(data.shape)
    np.fill_diagonal(data, 1.0)
    print(data)
    eigenvalues, eigenvectors = np.linalg.eig(data)
    print(eigenvalues.shape)
    print(eigenvectors.shape)
