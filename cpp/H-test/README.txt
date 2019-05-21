This file serves as a tool to test the implementation of operators and their action on states of the system. i.e.

Part of the Hamitlonian contains terms like 
$$\hat{H}_a = (\hat{a}_1 + \hat{a}_1^\dagger).$$
Such terms are implemented in the program, and during the execution of the main diagonalisation program, their action on all the basis vectors of the system is computed. So many operations like
$\bra{n_1, n_3, n_5; j_1, m_1, j_2, m_2} \hat{H}_a$
take place. The result of this operation is a superposition of many states, in case of $\hat{H}_a$ it is
$\bra{n_1, n_3, n_5; j_1, m_1, j_2, m_2} \hat{H}_a 
=
\bra{n_1+1, n_3, n_5; j_1, m_1, j_2, m_2}\sqrt{n_1+1}
+
\bra{n_1-1, n_3, n_5; j_1, m_1, j_2, m_2}\sqrt{n_1}.
$
This is what this program tests.
