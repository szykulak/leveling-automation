import sys

def gauss_elimination(matrix, number_of_unknowns): 
		solutions = [0 for i in
					 range(0, number_of_unknowns)]

		for i in range(number_of_unknowns):
			if matrix[i][i] == 0.0:
				sys.exit('Division by zero detected!')

			for j in range(i + 1, number_of_unknowns):
				ratio = matrix[j][i] / matrix[i][i]

				for k in range(number_of_unknowns + 1):
					matrix[j][k] = matrix[j][k] - ratio * matrix[i][k]

		solutions[number_of_unknowns - 1] = matrix[number_of_unknowns - 1][number_of_unknowns] / \
											matrix[number_of_unknowns - 1][number_of_unknowns - 1]

		for i in range(number_of_unknowns - 2, -1, -1):
			solutions[i] = matrix[i][number_of_unknowns]

			for j in range(i + 1, number_of_unknowns):
				solutions[i] = solutions[i] - matrix[i][j] * solutions[j]

			solutions[i] = solutions[i] / matrix[i][i]

		return solutions
