#####----- COMPOSITE MATERIAL PROPERTIES CALCULATOR WITH TKINTER -----#####
# By Francesco Sessa, Elaborato Costruzioni II

#----- SCRIPT SET UP -----#
# importing all the needed libraries
import os
import copy
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import colorsys
import csv
from datetime import datetime
from scipy.optimize import minimize

# Clean the terminal
os.system('cls')

#----- INTERACTIVE CALCULATOR FOR COMPOSITE LAMINATES -----#
class LaminateCalculator: 
    """Main class for the interactive calculator, a class is like a mould to create objects,
    but is also useful to organize the code in a better way"""

### Firstly the numerical and theoretical methods that the application follow to calculate, classificate and verify ###
    def calculate_everything(self, thickness_vector):
        """Calculate laminate properties using individual laminae with their own Q matrices"""
        if not self.laminae_list:
            tk.messagebox.showwarning("Warning", "No laminae defined! Please add laminae to the stack before calculating.")
            return
        try:
            # Read loads (GUI -> SI units)
            N_vec, M_vec, _ = self.get_loads_and_allowables()
            # Calculate Q matrices for each lamina
            Vector_Q_matrices = []
            rot2_Vector_Q_matrices = []
            rot3_Vector_Q_matrices = []
            rot4_Vector_Q_matrices = []
            rot5_Vector_Q_matrices = []

            z_positions = [0]  # Start from bottom
            
            # ps: We also calulate A2 and A3,A4,A5, A matrix calulated also in 2 different reference systems 
            for i, lamina in enumerate(self.laminae_list):
                # Convert units (MPa -> Pa)
                E_1 = lamina['E1'] * 1e6
                E_2 = lamina['E2'] * 1e6
                G_12 = lamina['G12'] * 1e6
                nu_12 = lamina['nu12']
                t_ply = thickness_vector[i]  # m
                angle = lamina['angle']
                angle2 = angle + 15
                angle3 = angle - 25
                angle4 = angle + 150
                angle5 = angle - 100
                
                # Calculate Q matrix for this lamina
                Q_11 = E_1/(1-(nu_12**2)*(E_2/E_1))
                Q_12 = (nu_12*E_2)/(1-(nu_12**2)*(E_2/E_1))
                Q_22 = E_2/(1-(nu_12**2)*(E_2/E_1))
                Q_66 = G_12
                Q = np.array([[Q_11, Q_12, 0],
                             [Q_12, Q_22, 0],
                             [0, 0, Q_66]])
                
                # Rotate Q matrix
                Theta_i = self.Theta(angle)
                Q_i = Theta_i @ Q @ Theta_i.T
                Vector_Q_matrices.append(Q_i)
                # Rotate Q matrix for other reference frames
                rot2_Theta_i = self.Theta(angle2)
                rot2_Q_i = rot2_Theta_i @ Q @ rot2_Theta_i.T
                rot2_Vector_Q_matrices.append(rot2_Q_i)
                rot3_Theta_i = self.Theta(angle3)
                rot3_Q_i = rot3_Theta_i @ Q @ rot3_Theta_i.T
                rot3_Vector_Q_matrices.append(rot3_Q_i)
                rot4_Theta_i = self.Theta(angle4)
                rot4_Q_i = rot4_Theta_i @ Q @ rot4_Theta_i.T
                rot4_Vector_Q_matrices.append(rot4_Q_i)
                rot5_Theta_i = self.Theta(angle5)
                rot5_Q_i = rot5_Theta_i @ Q @ rot5_Theta_i.T
                rot5_Vector_Q_matrices.append(rot5_Q_i)
                
                # Calculate z positions
                z_positions.append(z_positions[-1] + t_ply)
            
            # Center the z coordinates
            total_thickness = z_positions[-1]
            z = [z - total_thickness/2 for z in z_positions]
            
            # Calculate A, B, D matrices
            A = np.zeros((3,3))
            A2 = np.zeros((3,3))
            A3 = np.zeros((3,3))
            A4 = np.zeros((3,3))
            A5 = np.zeros((3,3))
            B = np.zeros((3,3))
            D = np.zeros((3,3))
            for i in range(len(self.laminae_list)):
                A = A + Vector_Q_matrices[i]*(z[i+1]-z[i])
                A2 = A2 + rot2_Vector_Q_matrices[i]*(z[i+1]-z[i])
                A3 = A3 + rot3_Vector_Q_matrices[i]*(z[i+1]-z[i])
                A4 = A4 + rot4_Vector_Q_matrices[i]*(z[i+1]-z[i])
                A5 = A5 + rot5_Vector_Q_matrices[i]*(z[i+1]-z[i])
                B = B + (1/2)*Vector_Q_matrices[i]*(z[i+1]**2 - z[i]**2)
                D = D + (1/3)*Vector_Q_matrices[i]*(z[i+1]**3 - z[i]**3)

            # Calculate inverse matrices a,b,c,d
            inv_A=np.linalg.inv(A)
            B_star = -inv_A @ B
            C_star = B @ inv_A
            D_star = D - B @ inv_A @ B
            inv_D_star = np.linalg.inv(D_star)
            a=inv_A - B_star @ inv_D_star @ C_star
            b=B_star @ inv_D_star
            c= b.T
            d=inv_D_star
            abcd=np.block([[a, b],
                [c, d]])
            # Calculate effective engineering constants of laminate
            Ex = 1/(a[0,0]*total_thickness)
            Ey = 1/(a[1,1]*total_thickness)
            Gxy = 1/(a[2,2]*total_thickness)
            nuxy = -a[0,1]/a[0,0]
            nuyx = -a[1,0]/a[1,1]
            Exb = 12/(d[0,0]*total_thickness**3)
            Eyb = 12/(d[1,1]*total_thickness**3)
            Dx = Exb*total_thickness**3/(12*(1 - nuxy*nuyx))
            Dy = Eyb*total_thickness**3/(12*(1 - nuxy*nuyx))
            #-- Calcolo Stress and Displacements (only if Load Condition is enabled) --#
            eps_i_global_list=[]; stress_i_global_list=[]
            eps_i_local_list=[];  stress_i_local_list=[]
            strain_FI_list=[];    stress_FI_list=[]
            if self.add_load_var.get():
                eps0_kappa = abcd @ np.concatenate((N_vec, M_vec))
                eps0=eps0_kappa[0:3]
                kappa=eps0_kappa[3:6]
            else:
                eps0 = np.zeros(3)
                kappa = np.zeros(3)
            # Get sequence of angles
            s1 = [lamina['angle'] for lamina in self.laminae_list]
            if self.add_load_var.get():
                for i in range(len(s1)):
                    # strain for each ply in global coordinates (mechanical part)
                    z_i = (z[i] + z[i+1])/2 
                    eps_i = eps0 + z_i * kappa
                    eps_i_global_list.append(eps_i)
                    # transformation matrices
                    T_i = self.Theta(s1[i])
                    T_i_inv = np.linalg.inv(T_i)
                    T_i_transp = T_i.T
                    # mechanical local strain
                    eps_i_local = T_i_transp @ eps_i
                    # thermal strain subtraction in local axes
                    lam = self.laminae_list[i]
                    a1 = lam.get('alpha_1', self.alpha1_micro_var.get()*1e-6)
                    a2 = lam.get('alpha_2', self.alpha2_micro_var.get()*1e-6)
                    dT = self.deltaT_var.get()
                    alpha_vec = np.array([a1, a2, 0.0])
                    eps_i_local_eff = eps_i_local - dT * alpha_vec
                    eps_i_local_list.append(eps_i_local)
                    # Compute stresses from local effective strain using Q (local)
                    # Local Q is the unrotated material matrix built earlier (Q)
                    # We need the local Q for this lamina; reconstruct from its properties
                    E_1 = lam['E1'] * 1e6
                    E_2 = lam['E2'] * 1e6
                    G_12 = lam['G12'] * 1e6
                    nu_12 = lam['nu12']
                    Q_11 = E_1/(1-(nu_12**2)*(E_2/E_1))
                    Q_12 = (nu_12*E_2)/(1-(nu_12**2)*(E_2/E_1))
                    Q_22 = E_2/(1-(nu_12**2)*(E_2/E_1))
                    Q_66 = G_12
                    Q_local = np.array([[Q_11, Q_12, 0],[Q_12, Q_22, 0],[0, 0, Q_66]])
                    stress_i_local = Q_local @ eps_i_local_eff
                    stress_i_local_list.append(stress_i_local)
                    # also provide global stress for completeness
                    stress_i = T_i @ stress_i_local
                    stress_i_global_list.append(stress_i)
                    # Per-lamina allowables (MPa -> Pa)
                    # Fallback to current defaults if preset-added lamina lacks allowables
                    s1lt = (lam.get('sigma_1_LT', self.sigma_1_LT_MPa.get()) or 0.0) * 1e6
                    s2lt = (lam.get('sigma_2_LT', self.sigma_2_LT_MPa.get()) or 0.0) * 1e6
                    s1lc = (lam.get('sigma_1_LC', self.sigma_1_LC_MPa.get()) or 0.0) * 1e6
                    s2lc = (lam.get('sigma_2_LC', self.sigma_2_LC_MPa.get()) or 0.0) * 1e6
                    tauL = (lam.get('tau_12_L', lam.get('tao_L', self.tao_L_MPa.get())) or 0.0) * 1e6
                    # Maximum Stress Criterion Failure Indexes 
                    stress_FI1_i = s1lt/stress_i_local[0] if stress_i_local[0] > 0 else s1lc/stress_i_local[0]
                    stress_FI2_i = s2lt/stress_i_local[1] if stress_i_local[1] > 0 else s2lc/stress_i_local[1]
                    stress_FI3_i = tauL/abs(stress_i_local[2]) if abs(stress_i_local[2]) > 1e-30 else float('inf')
                    stress_FI_list.append([stress_FI1_i, stress_FI2_i, stress_FI3_i])
                    # Per-lamina allowable strains using that lamina's moduli
                    E1i = (lam.get('E1', self.E1_var.get()) or 0.0) * 1e6
                    E2i = (lam.get('E2', self.E2_var.get()) or 0.0) * 1e6
                    G12i = (lam.get('G12', self.G12_var.get()) or 0.0) * 1e6
                    def sdiv(a,b):
                        return a/b if (b and abs(b)>1e-30) else float('inf')
                    eps1_lt = sdiv(s1lt, E1i)
                    eps2_lt = sdiv(s2lt, E2i)
                    eps1_lc = sdiv(s1lc, E1i)
                    eps2_lc = sdiv(s2lc, E2i)
                    gamma_l = sdiv(tauL, G12i)
                    epsi_FI1_i = eps1_lt/eps_i_local[0] if eps_i_local[0] > 0 else eps1_lc/eps_i_local[0]
                    epsi_FI2_i = eps2_lt/eps_i_local[1] if eps_i_local[1] > 0 else eps2_lc/eps_i_local[1]
                    epsi_FI3_i = gamma_l/abs(eps_i_local[2]) if abs(eps_i_local[2]) > 1e-30 else float('inf')
                    strain_FI_list.append([epsi_FI1_i, epsi_FI2_i, epsi_FI3_i])
            
            results={
                'Vector_Q_matrices': Vector_Q_matrices,
                'A': A,
                'B': B,
                'D': D,
                's1': s1,
                'A2': A2,
                'A3': A3,
                'A4': A4,
                'A5': A5,
                'Ex': Ex,
                'Ey': Ey,
                'Gxy': Gxy,
                'nuxy': nuxy,
                'nuyx': nuyx,
                'Exb': Exb,
                'Eyb': Eyb,
                'Dx': Dx,
                'Dy': Dy,
                't': total_thickness,
                'z': z,
                'eps0': eps0,
                'kappa': kappa,
                'eps_i_global_list': eps_i_global_list,
                'stress_i_global_list': stress_i_global_list,
                'eps_i_local_list': eps_i_local_list,
                'stress_i_local_list': stress_i_local_list,
                'strain_FI_list': strain_FI_list,
                'stress_FI_list': stress_FI_list,
            }
            return results
        except Exception as e:
            tk.messagebox.showerror("Error", f"Error in calculations: {str(e)}")
            return None

    def Theta(self, angle): # Method to calculate the rotation matrix
        """Calculate the transformation matrix for a given angle in degrees"""
        m = np.cos(np.radians(angle))
        n = np.sin(np.radians(angle))
        return np.array([[m**2, n**2, -2*m*n],
                        [n**2, m**2, 2*m*n],
                        [m*n, -m*n, m**2 - n**2]])
    
    def fmt_small(self, x, precision=3, zero="0"):
        """Format numbers: show zero string if |x| < 1e-7 else scientific with given precision."""
        try:
            xv = float(x)
        except Exception:
            return str(x)
        if abs(xv) < 1e-7:
            return zero
        fmt = "{:." + str(precision) + "e}"
        return fmt.format(xv)
    
    def classify_laminate(self, angles, A_matrix, B_matrix, D_matrix, A2, A3, A4,A5):
        """Classify the laminate based on its special properties"""
        # Normalize angles to be within -180 to 180 degrees
        angle_counts = -1
        for angle in angles:
           angle_counts += 1
           if angle>180:
                angles[angle_counts] = angle - 360

        special_cases = []
        is_symmetric = self.is_symmetric(angles)
        is_antisymmetric = self.is_antisymmetric(angles)
        is_crossply = self.is_crossply(angles)
        is_balanced = self.is_balanced(angles)
        is_quasi_isotropic = self.is_quasi_isotropic(A_matrix,A2,A3,A4,A5)
        # 1. Check if symmetric and if symmetric cross-ply and if antisymmetric and if antisymmetric cross-ply
        if is_symmetric:
            if is_crossply:
                special_cases.append("Symmetric Cross-Ply")
            else:
                special_cases.append("Symmetric")
        elif is_crossply and self.is_zero_ninety(angles):
            special_cases.append("Antisymmetric Cross-Ply")
        elif (is_antisymmetric):
            if (is_crossply) :
                special_cases.append("Antisymmetric Cross-Ply")
            else:
                special_cases.append("Antisymmetric")
        elif (is_crossply):
            special_cases.append("Cross-Ply")
        # 2. Check if balanced
        if is_balanced:
            special_cases.append("Balanced")
        # 3. Check if quasi-isotropic
        if is_quasi_isotropic:
            special_cases.append("Quasi-Isotropic")
        return special_cases
    
    def is_symmetric(self, angles):
        """Check if laminate is symmetric"""
        n = len(angles)
        for i in range(n//2):
            lamina_i = self.laminae_list[0]
            lamina_alt= self.laminae_list[n-1-i]  
            if angles[i]==-90 or angles[n-1-i]==-90:
                angles[i]=90
                angles[n-1-i]=90
            if angles[i] != angles[n-1-i]:
                return False
            if (lamina_i['E1'] != lamina_alt['E1'] or
                lamina_i['E2'] != lamina_alt['E2'] or
                lamina_i['G12'] != lamina_alt['G12'] or
                lamina_i['nu12'] != lamina_alt['nu12'] or
                lamina_i['thickness'] != lamina_alt['thickness']):
                return False
        return True
    
    def is_antisymmetric(self, angles):
        """Check if laminate is antisymmetric"""
        n = len(angles)
        for i in range(n//2):
            lamina_i = self.laminae_list[i]
            lamina_alt= self.laminae_list[n-1-i]
            if (lamina_i['E1'] != lamina_alt['E1'] or
                lamina_i['E2'] != lamina_alt['E2'] or
                lamina_i['G12'] != lamina_alt['G12'] or
                lamina_i['nu12'] != lamina_alt['nu12'] or
                lamina_i['thickness'] != lamina_alt['thickness']):
                return False
            if abs(angles[i]) != 90 and angles[i] != 0:
                if angles[i] != -angles[n-1-i]:
                    return False
            elif abs(angles[i]) == 90:
                if abs(angles[n-1-i]) != 90:
                    return False
            elif angles[i] == 0:
                if angles[n-1-i] != 0:
                    return False
        return True
    
    def is_balanced(self, angles):
        """Check if laminate is balanced"""
        positive_angles = [angle for angle in angles if angle > 0 and abs(angle) != 90]
        negative_angles = [angle for angle in angles if angle < 0 and abs(angle) != 90]
        positive_angles.sort()
        negative_angles.sort(reverse=True)
        if len(positive_angles) != len(negative_angles):
            return False
        for i in range(len(positive_angles)):
            if positive_angles[i] + negative_angles[i]:
                return False
        ang_90_count = angles.count(90) + angles.count(-90)
        ang_0_count = angles.count(0)
        if ang_90_count % 2 != 0 or ang_0_count % 2 != 0:
            return False
        return True
    
    def is_crossply(self,angles):
        """Check if cross-ply"""
        if not self.is_same_lamine_property():
            return False
        for angle in angles:
            if angle not in [0, 90, -90]:
                return False
        return True
    
    def is_zero_ninety(self, angles):
        """Check if antisymmetric cross-ply, so made of sequence of [0/90] degrees only
        remember to check before if it is cross-ply"""
        n = len(angles)
        first_angle = angles[0]
        for i in range(n):
            if not i%2:
                if first_angle == 0:
                    if angles[i]==0 and angles[i+1] == 90:
                        continue
                    else:
                        return False
                elif first_angle == 90:
                    if angles[i] == 90:
                        if angles[i+1] == 0 and angles[i+1] == 0:
                            continue
                        else:
                            return False
        return True
               
    def is_quasi_isotropic(self, A_matrix, A2, A3, A4, A5):
        """Check if laminate is quasi-isotropic"""
        tolerance = 1e-4
        if not self.is_same_lamine_property():
            return False
        for i in range(3):
            for j in range(3):
                if abs(A_matrix[i,j] - A2[i,j]) > tolerance or \
                   abs(A_matrix[i,j] - A3[i,j]) > tolerance or \
                   abs(A_matrix[i,j] - A4[i,j]) > tolerance or \
                   abs(A_matrix[i,j] - A5[i,j]) > tolerance:
                    return False
        return True

    def is_same_lamine_property(self):
        """Check if all laminae have the same material properties and thickness"""
        if len(self.laminae_list) > 1:
            first_lamina = self.laminae_list[0]
            reference_E1 = first_lamina['E1']
            reference_E2 = first_lamina['E2']
            reference_G12 = first_lamina['G12']
            reference_nu12 = first_lamina['nu12']
            reference_thickness = first_lamina['thickness']
            # Check if all laminae have the same material properties and thickness
            for lamina in self.laminae_list[1:]:
                if (lamina['E1'] != reference_E1 or
                    lamina['E2'] != reference_E2 or
                    lamina['G12'] != reference_G12 or
                    lamina['nu12'] != reference_nu12 or
                    lamina['thickness'] != reference_thickness):
                    return False
        return True

    def assess_failure(self, results):
        """Assess failure; return dict with pass flag, min SM stress/strain and failed plies info."""
        s1 = results['s1']
        stress_FI_list = results['stress_FI_list']
        strain_FI_list = results['strain_FI_list']
        failures = []  # list of dicts: {index, reasons: [(type, dir, FI)]}
        min_sm_stress = float('inf')
        min_sm_strain = float('inf')
        for i in range(len(s1)):
            ply_fail = []
            # stress criteria
            min_stress_FI = min(stress_FI_list[i])
            min_sm_stress = min(min_sm_stress, min_stress_FI - 1)
            for j, comp in enumerate(['1','2','12']):
                if stress_FI_list[i][j] < 1:
                    ply_fail.append(('stress', comp, float(stress_FI_list[i][j])))
            # strain criteria
            min_strain_FI = min(strain_FI_list[i])
            min_sm_strain = min(min_sm_strain, min_strain_FI - 1)
            for j, comp in enumerate(['1','2','12']):
                if strain_FI_list[i][j] < 1:
                    ply_fail.append(('strain', comp, float(strain_FI_list[i][j])))
            if ply_fail:
                failures.append({'index': i, 'reasons': ply_fail})
        passed = len(failures) == 0
        return {
            'passed': passed,
            'min_sm_stress': float(min_sm_stress if min_sm_stress != float('inf') else 0.0),
            'min_sm_strain': float(min_sm_strain if min_sm_strain != float('inf') else 0.0),
            'failures': failures
        }

    def optimize_laminate(self):
        """Optimize the laminate thicknesses to meet the total thickness and safety margin requirements"""
        # Snapshot the current laminate so we can restore the exact ply set when applying optimized thicknesses
        try:
            self._laminate_snapshot = copy.deepcopy(self.laminae_list)
        except Exception:
            self._laminate_snapshot = None
        # Build constraints and bounds from GUI
       
        total_thickness_target_m = self.total_thickness_target_mm_var.get()/1000.0
        lam_boundary_thickness = (self.lam_min_thick_mm_var.get()/1000.0,
                                    self.lam_max_thick_mm_var.get()/1000.0)
        bounds_list = [ lam_boundary_thickness for _ in range(len(self.thickness_vector))]
        # Build an initial guess closer to the target total thickness to help convergence
        nplies = max(1, len(self.thickness_vector))
        uniform_guess = total_thickness_target_m / nplies
        uniform_guess = min(max(uniform_guess, lam_boundary_thickness[0]), lam_boundary_thickness[1])
        t0_vector= [uniform_guess for _ in range (len(self.thickness_vector))]
        # if there is a failure
        if self.failed_plies:
            constraints_list=[
            {'type': 'ineq', 'fun': lambda x: total_thickness_target_m-sum(x)},  # vincolo di spessore totale
            {'type': 'ineq', 'fun': self.minim_SM_constraint} # vincolo di sicurezza minima
            ]
            if self.same_thickness_var.get():
                constraints_list.append({'type': 'eq', 'fun': self.same_thickness_constraint})
            fun_criteria=self.obj_stress_criteria if self.stress_criteria_var.get() else self.obj_strain_criteria
            optimiz_result = minimize(
                fun=fun_criteria, 
                x0=t0_vector, 
                bounds=bounds_list,
                constraints=constraints_list,
                method='SLSQP',
                options={'maxiter': 500, 'ftol': 1e-9, 'disp': False}
            )
        else:
            constraints_list=[
            {'type': 'ineq', 'fun': lambda x: total_thickness_target_m-sum(x)},  # vincolo di spessore totale
             ]
            if self.same_thickness_var.get():
                constraints_list.append({'type': 'eq', 'fun': self.same_thickness_constraint})
           
            if self.stress_criteria_var.get():
                constraints_list.append({'type': 'ineq', 'fun': lambda x: -self.obj_stress_criteria(x)}) # vincolo di sicurezza minima
            else:
                constraints_list.append({'type': 'ineq', 'fun': lambda x: -self.obj_strain_criteria(x)}) # vincolo di sicurezza minima
        optimiz_result = minimize(
            fun=self.obj_decrease_SM, 
            x0=t0_vector, 
            bounds=bounds_list,
            constraints=constraints_list,
            method='SLSQP',
            options={'maxiter': 500, 'ftol': 1e-9, 'disp': False}
        )

        # Save summary for the Optimization Results tab
        self.optimization_result = {
            'success': bool(optimiz_result.success),
            'message': str(optimiz_result.message),
            'nit': int(optimiz_result.nit),
            'x': np.array(optimiz_result.x).tolist(),
            'fun': float(optimiz_result.fun),
        }
        new_thick_vector=optimiz_result.x
        new_result=self.calculate_everything(new_thick_vector)
        if self.stress_criteria_var.get():
            min_SM=min([min(new_result['stress_FI_list'][i]) for i in range(len(new_result['stress_FI_list']))])-1
        else:
            min_SM=min([min(new_result['strain_FI_list'][i]) for i in range(len(new_result['strain_FI_list']))])-1
        new_total_thick=float(sum(new_thick_vector))
        self.optimization_result['min_SM_after'] = float(min_SM)
        self.optimization_result['total_thickness_after'] = new_total_thick
        # Update tab
        self.update_optimization_tab()

    def minim_SM_constraint(self,thickness_vector):
        """Constraint function to ensure minimum safety margin is met during optimization"""
        if self.stress_criteria_var.get():
            stress_FI_list=self.calculate_everything(thickness_vector)['stress_FI_list']
            min_SM=min([min(stress_FI_list[i]) for i in range(len(stress_FI_list))])-1
        else:
            strain_FI_list=self.calculate_everything(thickness_vector)['strain_FI_list']
            min_SM=min([min(strain_FI_list[i]) for i in range(len(strain_FI_list))])-1
        return self.SM_wanted_var.get()-min_SM
    
    def same_thickness_constraint(self, t_vec):
        """Constraint function to ensure all laminae have the same thickness"""
        return [t_vec[i] - t_vec[0] for i in range(1, len(t_vec))]
    
    def obj_stress_criteria(self,t_vec):
        new_stress_FI_list=self.calculate_everything(t_vec)['stress_FI_list']
        # objective value = somma di tutti gli spessori delle lamine che non rispettano il Safety Margin richiesto
        # add small epsilon in denominators to avoid infinities/steep gradients
        eps = 1e-12
        objective_value = sum([
            t_vec[i]/max(min(new_stress_FI_list[i]), eps)
            for i in range(len(t_vec))
            if min(new_stress_FI_list[i])<(self.SM_wanted_var.get()+1)
        ])
        return objective_value

    def obj_strain_criteria(self,t_vec):
        new_strain_FI_list=self.calculate_everything(t_vec)['strain_FI_list']
        # objective value = somma di tutti gli spessori delle lamine che non rispettano il Safety Margin richiesto
        eps = 1e-12
        objective_value = sum([
            t_vec[i]/max(min(new_strain_FI_list[i]), eps)
            for i in range(len(t_vec))
            if min(new_strain_FI_list[i])<(self.SM_wanted_var.get()+1)
        ])
        return objective_value
    
    def obj_decrease_SM(self,t_vec):
        # Use smooth squared penalty instead of absolute value to aid SLSQP line search
        v = self.minim_SM_constraint(t_vec)
        return v*v
    
    def apply_optimized_thicknesses(self):
        """Apply optimized thickness vector to current laminate as input (convert m -> mm)."""
        if not self.optimization_result or 'x' not in self.optimization_result:
            tk.messagebox.showinfo("Optimization", "No optimized thicknesses available. Run optimization first.")
            return
        x_m = np.array(self.optimization_result['x'])
        if len(self.laminae_list) == 0 and (getattr(self, '_laminate_snapshot', None) is None):
            tk.messagebox.showwarning("Laminate", "No laminae defined to apply thicknesses.")
            return
        # Prefer restoring the exact snapshot of the laminate (plies, order, angles, materials)
        snap = getattr(self, '_laminate_snapshot', None)
        if snap is not None and len(snap) == len(x_m):
            self.laminae_list = copy.deepcopy(snap)
            for i in range(len(self.laminae_list)):
                self.laminae_list[i]['thickness'] = float(x_m[i]*1000.0)
        else:
            # Fallback: update existing list up to the available length
            n = min(len(self.laminae_list), len(x_m))
            for i in range(n):
                self.laminae_list[i]['thickness'] = float(x_m[i]*1000.0)  # store in mm
        # Update thickness vector in meters and refresh UI
        self.thickness_vector = [lam['thickness']/1000.0 for lam in self.laminae_list]
        self.update_laminae_display()
        # Optionally recalc results with new thicknesses
        self.update_results_display()
    
### Constructor of the class that start the call of the other methods and define the properties ###
    def __init__(self): 
        """ Constructor of the class, initializes the GUI,
        here we define the properties: root and input variables"""

        # property root (main window)
        self.root = tk.Tk() # create the main window
        self.root.title("ICLC - Francesco Sessa")
        # Interactive Composite Laminate Calculator=ICLC
        self.root.geometry("1400x900")
        self.root.configure(bg='lightgray')

        # Other Properties of the main window are the Tkinter variables for inputs
        # Moduli now requested in MPa (previous defaults were in GPa); convert defaults GPa->MPa
        self.E1_var = tk.DoubleVar(value=125000.0) # Young modulus in MPa
        self.E2_var = tk.DoubleVar(value=12500.0)  # MPa
        self.G12_var = tk.DoubleVar(value=6890.0)  # MPa (6.89 GPa)
        self.nu12_var = tk.DoubleVar(value=0.38)
        self.angle_var = tk.DoubleVar(value=0.0)
        self.thickness_var = tk.DoubleVar(value=0.2) # thickness of the single layer in mm
        self.laminae_list = []  # List to store individual laminae data
        self.thickness_vector = []  # List to store thicknesses of laminae
        self.result = {}  # Dictionary to store calculation results
        # Thermal expansion coefficients input as micro-scale (user enters -1, 50 meaning -1e-6, 50e-6)
        self.alpha1_micro_var = tk.DoubleVar(value=-1.0)
        self.alpha2_micro_var = tk.DoubleVar(value=50.0)
        # --- New GUI Variables ---
        # Load conditions (shown when checkbox enabled)
        self.add_load_var = tk.BooleanVar(value=False)
        self.Nx_var = tk.DoubleVar(value=8500.0)
        self.Ny_var = tk.DoubleVar(value=17000.0)
        self.Nxy_var = tk.DoubleVar(value=16000.0)
        self.Mx_var = tk.DoubleVar(value=100.0)
        self.My_var = tk.DoubleVar(value=0.0)
        self.Mxy_var = tk.DoubleVar(value=0.0)
        # Temperature change (°C) for thermal strain (only applied if Load Condition enabled)
        self.deltaT_var = tk.DoubleVar(value=0.0)
        # Allowables (MPa in GUI, converted to Pa in code)
        self.sigma_1_LT_MPa = tk.DoubleVar(value=434.0)
        self.sigma_2_LT_MPa = tk.DoubleVar(value=41.4)
        self.sigma_1_LC_MPa = tk.DoubleVar(value=-331.0)
        self.sigma_2_LC_MPa = tk.DoubleVar(value=-27.6)
        self.tao_L_MPa      = tk.DoubleVar(value=34.5)
        # Optimization options
        self.optimize_var = tk.BooleanVar(value=False)
        self.SM_wanted_var = tk.DoubleVar(value=0.2)
        self.same_thickness_var = tk.BooleanVar(value=True)
        self.total_thickness_target_mm_var = tk.DoubleVar(value=4.6)  # mm
        self.lam_min_thick_mm_var = tk.DoubleVar(value=0.1)
        self.lam_max_thick_mm_var = tk.DoubleVar(value=0.5)
        # Criteria toggles (mutually exclusive; default to Stress)
        self.stress_criteria_var = tk.BooleanVar(value=True)
        self.strain_criteria_var = tk.BooleanVar(value=False)
        # Failure info
        self.failed_plies = {}
        # Optimization result store
        self.optimization_result = None
        # Clipboard for copy/cut/paste of laminae
        self._lamina_clipboard = []

        self.setup_interface() # Inside the constructor we call the method to set up the interface
        # A method is like a function that belongs to a class called with the syntax NameofC self.method_name()
        self.load_preset_1() # Load default preset [0/45/-45/90]s

### From here Methods for setting up the interface, managing user inputs and outputs ###
    def setup_interface(self):
        """ In this method we set up the interface of the calculator, we divide the window in two columns,
        to do that we create a main frame that holds two frames: left for inputs, right for results tabs"""

        # Main Frame
        main_frame = tk.Frame(self.root, bg='lightgray')
        main_frame.pack(fill='both', expand=True, padx=10, pady=10) # pack the main frame to fill the root window
        # pack is a geometry manager that organizes widgets in blocks before placing them in the parent widget
        # after creating a tkinter widget, you need to call a geometry manager method to display it in the window

        # Left column - Input parameters (scrollable, fully contained)
        input_frame = tk.LabelFrame(main_frame, text="Laminate Design",
                                   font=("Arial", 14, "bold"), bg='white', pady=0, width=460)
        # Fix the input panel width so inner widgets don't force it to overflow
        input_frame.pack_propagate(False)
        input_frame.pack(side='left', fill='y', padx=(0, 10))
        input_canvas = tk.Canvas(input_frame, bg='white', highlightthickness=0)
        input_scrollbar = ttk.Scrollbar(input_frame, orient="vertical", command=input_canvas.yview)
        self.input_content = tk.Frame(input_canvas, bg='white')
        # Slightly narrower so everything stays within the input panel
        input_canvas.configure(yscrollcommand=input_scrollbar.set, width=420)
        input_canvas.pack(side="left", fill="both", expand=False)
        input_scrollbar.pack(side="right", fill="y")
        # Keep the inner window id to dynamically clamp width on resize
        self._input_window_id = input_canvas.create_window((0, 0), window=self.input_content, anchor="nw")
        def _on_input_configure(event):
            # Update scrollregion whenever the content changes size
            input_canvas.configure(scrollregion=input_canvas.bbox("all"))
        self.input_content.bind("<Configure>", _on_input_configure)
        def _on_frame_resize(event):
            # Constrain the inner content width to the available canvas width minus scrollbar
            avail = max(100, event.width - 16)
            input_canvas.configure(width=avail)
            try:
                input_canvas.itemconfigure(self._input_window_id, width=avail)
            except Exception:
                pass
        input_frame.bind("<Configure>", _on_frame_resize)

        # Initialize input section
        self.create_input_section(self.input_content) # method to initialize the input section, 
        # inside it all the widgets for inputs are created

        # Right column - Results
        results_frame = tk.Frame(main_frame, bg='lightgray')
        results_frame.pack(side='right', fill='both', expand=True)
        
        # The Right Column is divided in more tabs,
        # we need to create a notebook to hold results tabs
        self.notebook = ttk.Notebook(results_frame)
        self.notebook.pack(fill='both', expand=True)
        # Tab for Laminate Visual Representation
        self.visual_frame = tk.Frame(self.notebook, bg='white')
        self.notebook.add(self.visual_frame, text="Laminate Visual Representation")
        # Tab for the matrices
        self.matrices_frame = tk.Frame(self.notebook, bg='white')
        self.notebook.add(self.matrices_frame, text="Matrices Q, A, B, D")
        # Tab for the properties (moved before Stress & Strain)
        self.properties_frame = tk.Frame(self.notebook, bg='white')
        self.notebook.add(self.properties_frame, text="Engineering Properties")
        # Tab for Stress & Strain (now after properties and before optimization)
        self.checks_frame = tk.Frame(self.notebook, bg='white')
        self.notebook.add(self.checks_frame, text="Stress & Strain")
        # Tab for Optimization Results
        self.optim_frame = tk.Frame(self.notebook, bg='white')
        self.notebook.add(self.optim_frame, text="Optimization Results")
        # Tab for Credits (last)
        self.credits_frame = tk.Frame(self.notebook, bg='white')
        self.notebook.add(self.credits_frame, text="Credits")
        self._setup_credits_tab()
        
        # Initialize result widgets
        self.setup_results_widgets() # method to initialize the results section, 
        # inside it all the widgets for results are created
        # Apply numeric validation to all input entries (commas/non-numeric)
        try:
            self._attach_numeric_validation_all()
        except Exception:
            pass

### Managing the input section and the inputs ###
    def create_input_section(self, parent):
        """ Here we define all the widgets for the input section in the left column of the GUI"""
        # Removed header label to save vertical space
        add_frame = tk.LabelFrame(parent, text="Add New Lamina", font=("Arial", 12, "bold"), bg='white')
        add_frame.pack(fill='x', padx=5, pady=5)
        props_frame = tk.Frame(add_frame, bg='white')
        props_frame.pack(fill='x', padx=5, pady=5)
        row1_frame = tk.Frame(props_frame, bg='white'); row1_frame.pack(fill='x', pady=2)
        tk.Label(row1_frame, text="E₁ (MPa):", font=("Arial", 10), bg='white').pack(side='left')
        tk.Entry(row1_frame, textvariable=self.E1_var, font=("Arial", 10), width=8).pack(side='left', padx=4)
        tk.Label(row1_frame, text="E₂ (MPa):", font=("Arial", 10), bg='white').pack(side='left', padx=(8,0))
        tk.Entry(row1_frame, textvariable=self.E2_var, font=("Arial", 10), width=8).pack(side='left', padx=4)
        tk.Label(row1_frame, text="G₁₂ (MPa):", font=("Arial", 10), bg='white').pack(side='left', padx=(8,0))
        tk.Entry(row1_frame, textvariable=self.G12_var, font=("Arial", 10), width=8).pack(side='left', padx=4)
        row2_frame = tk.Frame(props_frame, bg='white'); row2_frame.pack(fill='x', pady=2)
        tk.Label(row2_frame, text="ν₁₂:", font=("Arial", 10), bg='white').pack(side='left')
        tk.Entry(row2_frame, textvariable=self.nu12_var, font=("Arial", 10), width=8).pack(side='left', padx=4)
        tk.Label(row2_frame, text="t (mm):", font=("Arial", 10), bg='white').pack(side='left', padx=(8,0))
        tk.Entry(row2_frame, textvariable=self.thickness_var, font=("Arial", 10), width=8).pack(side='left', padx=4)
        tk.Label(row2_frame, text="Angle (°):", font=("Arial", 10), bg='white').pack(side='left', padx=(8,0))
        tk.Entry(row2_frame, textvariable=self.angle_var, font=("Arial", 10), width=8).pack(side='left', padx=4)
        # Row 2b: thermal expansion coefficients (micro-scale user inputs; stored converted later)
        row2b_frame = tk.Frame(props_frame, bg='white'); row2b_frame.pack(fill='x', pady=2)
        tk.Label(row2b_frame, text="α_1 (1E-6):", font=("Arial", 10), bg='white').pack(side='left')
        tk.Entry(row2b_frame, textvariable=self.alpha1_micro_var, font=("Arial", 10), width=6).pack(side='left', padx=4)
        tk.Label(row2b_frame, text="α_2 (1E-6):", font=("Arial", 10), bg='white').pack(side='left', padx=(8,0))
        tk.Entry(row2b_frame, textvariable=self.alpha2_micro_var, font=("Arial", 10), width=6).pack(side='left', padx=4)
        allow_hdr = tk.Label(props_frame, text="Allowables (MPa)", font=("Arial", 10, 'bold'), bg='white')
        allow_hdr.pack(anchor='w', padx=5, pady=(2,0))
        row3a = tk.Frame(props_frame, bg='white'); row3a.pack(fill='x', pady=(1,0))
        tk.Label(row3a, text="σ1_LT:", font=("Arial", 10), bg='white').pack(side='left')
        tk.Entry(row3a, textvariable=self.sigma_1_LT_MPa, font=("Arial", 10), width=6).pack(side='left', padx=2)
        tk.Label(row3a, text="σ2_LT:", font=("Arial", 10), bg='white').pack(side='left', padx=(6,0))
        tk.Entry(row3a, textvariable=self.sigma_2_LT_MPa, font=("Arial", 10), width=6).pack(side='left', padx=2)
        tk.Label(row3a, text="τ12_L:", font=("Arial", 10), bg='white').pack(side='left', padx=(6,0))
        tk.Entry(row3a, textvariable=self.tao_L_MPa, font=("Arial", 10), width=6).pack(side='left', padx=2)
        row3b = tk.Frame(props_frame, bg='white'); row3b.pack(fill='x', pady=(1,0))
        tk.Label(row3b, text="σ1_LC:", font=("Arial", 10), bg='white').pack(side='left')
        tk.Entry(row3b, textvariable=self.sigma_1_LC_MPa, font=("Arial", 10), width=6).pack(side='left', padx=2)
        tk.Label(row3b, text="σ2_LC:", font=("Arial", 10), bg='white').pack(side='left', padx=(6,0))
        tk.Entry(row3b, textvariable=self.sigma_2_LC_MPa, font=("Arial", 10), width=6).pack(side='left', padx=2)
        # Button to add lamina
        add_button = tk.Button(add_frame, text="Add Lamina", command=self.add_lamina,
                              font=("Arial", 10, "bold"), bg='green', fg='white')
        add_button.pack(pady=10)
        # Frame for displaying laminae list
        list_frame = tk.LabelFrame(parent, text="Current Laminate Stack", 
                                  font=("Arial", 12, "bold"), bg='white')
        # Keep this frame within the input panel width (no horizontal growth)
        list_frame.pack(fill='x', expand=False, padx=5, pady=5)
        # Container just for the laminate stack table + its scrollbars
        tree_container = tk.Frame(list_frame, bg='white', height=190)
        tree_container.pack(fill='x', expand=False, padx=5, pady=(5,0))
        tree_container.grid_propagate(False)  # keep a compact, fixed-height area
        tree_container.grid_rowconfigure(0, weight=1)
        tree_container.grid_columnconfigure(0, weight=1)
        # Create treeview for laminae list
        columns = ('layer', 'E1', 'E2', 'G12', 'nu12', 'thickness', 'angle', 'a1', 'a2', 's1LT', 's2LT', 's1LC', 's2LC', 't12L')
        self.laminae_tree = ttk.Treeview(tree_container, columns=columns, show='headings', height=8)
        # Configure columns
        self.laminae_tree.heading('layer', text='Layer')
        self.laminae_tree.heading('E1', text='E₁ (MPa)')
        self.laminae_tree.heading('E2', text='E₂ (MPa)')
        self.laminae_tree.heading('G12', text='G₁₂ (MPa)')
        self.laminae_tree.heading('nu12', text='ν₁₂')
        self.laminae_tree.heading('thickness', text='t (mm)')
        self.laminae_tree.heading('angle', text='θ (°)')
        self.laminae_tree.heading('a1', text='α1 (×10⁻⁶)')
        self.laminae_tree.heading('a2', text='α2 (×10⁻⁶)')
        self.laminae_tree.heading('s1LT', text='σ1_LT (MPa)')
        self.laminae_tree.heading('s2LT', text='σ2_LT (MPa)')
        self.laminae_tree.heading('s1LC', text='σ1_LC (MPa)')
        self.laminae_tree.heading('s2LC', text='σ2_LC (MPa)')
        self.laminae_tree.heading('t12L', text='τ12_L (MPa)')
        self.laminae_tree.column('layer', width=48, anchor='center', stretch=False)
        self.laminae_tree.column('E1', width=64, anchor='center', stretch=False)
        self.laminae_tree.column('E2', width=64, anchor='center', stretch=False)
        self.laminae_tree.column('G12', width=64, anchor='center', stretch=False)
        self.laminae_tree.column('nu12', width=58, anchor='center', stretch=False)
        self.laminae_tree.column('thickness', width=58, anchor='center', stretch=False)
        self.laminae_tree.column('angle', width=56, anchor='center', stretch=False)
        self.laminae_tree.column('a1', width=70, anchor='center', stretch=False)
        self.laminae_tree.column('a2', width=70, anchor='center', stretch=False)
        self.laminae_tree.column('s1LT', width=86, anchor='center', stretch=False)
        self.laminae_tree.column('s2LT', width=86, anchor='center', stretch=False)
        self.laminae_tree.column('s1LC', width=86, anchor='center', stretch=False)
        self.laminae_tree.column('s2LC', width=86, anchor='center', stretch=False)
        self.laminae_tree.column('t12L', width=86, anchor='center', stretch=False)
        # Add dedicated scrollbars so this wide widget never gets cut off
        xscroll = ttk.Scrollbar(tree_container, orient='horizontal', command=self.laminae_tree.xview)
        yscroll = ttk.Scrollbar(tree_container, orient='vertical', command=self.laminae_tree.yview)
        self.laminae_tree.configure(xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
        # Layout with grid to keep a small side vertical scrollbar only for the table
        self.laminae_tree.grid(row=0, column=0, sticky='nsew')
        yscroll.grid(row=0, column=1, sticky='ns')
        xscroll.grid(row=1, column=0, columnspan=2, sticky='ew')
        # Buttons for managing laminae
        buttons_frame = tk.Frame(list_frame, bg='white')
        buttons_frame.pack(fill='x', padx=5, pady=5)
        remove_btn = tk.Button(buttons_frame, text="Remove Selected", command=self.remove_lamina,
                              font=("Arial", 9), bg='red', fg='white', width=12)
        remove_btn.pack(side='left', padx=2)
        clear_btn = tk.Button(buttons_frame, text="Clear All", command=self.clear_laminae,
                             font=("Arial", 9), bg='orange', fg='white', width=12)
        clear_btn.pack(side='left', padx=2)
        move_up_btn = tk.Button(buttons_frame, text="Move Up", command=self.move_lamina_up,
                               font=("Arial", 9), bg='blue', fg='white', width=12)
        move_up_btn.pack(side='left', padx=2)
        move_down_btn = tk.Button(buttons_frame, text="Move Down", command=self.move_lamina_down,
                                 font=("Arial", 9), bg='blue', fg='white', width=12)
        move_down_btn.pack(side='left', padx=2)
        # Clipboard controls (Copy / Cut / Paste) just below
        clipboard_frame = tk.Frame(list_frame, bg='white')
        clipboard_frame.pack(fill='x', padx=5, pady=(0,6))
        copy_btn = tk.Button(clipboard_frame, text="Copy", command=self.copy_laminae,
                             font=("Arial", 9), bg='#cccccc', width=10)
        copy_btn.pack(side='left', padx=2)
        cut_btn = tk.Button(clipboard_frame, text="Cut", command=self.cut_laminae,
                            font=("Arial", 9), bg='#bbbbbb', width=10)
        cut_btn.pack(side='left', padx=2)
        paste_btn = tk.Button(clipboard_frame, text="Paste", command=self.paste_laminae,
                              font=("Arial", 9), bg='#a6e3a1', width=10)
        paste_btn.pack(side='left', padx=2)
    # (moved) Main calculation buttons will be placed at the bottom of the input section
        # Quick preset buttons
        preset_frame = tk.LabelFrame(parent, text="Quick Presets", font=("Arial", 10, "bold"), bg='white')
        preset_frame.pack(fill='x', padx=5, pady=5)
        preset_buttons_frame = tk.Frame(preset_frame, bg='white')
        preset_buttons_frame.pack(fill='x', padx=5, pady=5)

        preset1_btn = tk.Button(preset_buttons_frame, text="[0/45/-45/90]s", 
                               command=self.load_preset_1, font=("Arial", 9), bg='lightblue', width=12)
        preset1_btn.pack(side='left', padx=2)

        preset2_btn = tk.Button(preset_buttons_frame, text="[0/90]s", 
                               command=self.load_preset_2, font=("Arial", 9), bg='lightblue', width=12)
        preset2_btn.pack(side='left', padx=2)
        
        preset3_btn = tk.Button(preset_buttons_frame, text="[±45]s", 
                               command=self.load_preset_3, font=("Arial", 9), bg='lightblue', width=12)
        preset3_btn.pack(side='left', padx=2)

        # --- Load condition & Allowables section ---
        load_chk_frame = tk.Frame(parent, bg='white')
        load_chk_frame.pack(fill='x', padx=5, pady=(10,0))
        load_chk = ttk.Checkbutton(load_chk_frame, text="Add Load Condition", variable=self.add_load_var, command=self.toggle_load_section)
        load_chk.pack(anchor='w')
        self.load_section = tk.LabelFrame(parent, text="Loads (N/m)", font=("Arial", 11, "bold"), bg='white')
        # inside load_section (initially hidden until checkbox)
        # Loads
        loads_row1 = tk.Frame(self.load_section, bg='white'); loads_row1.pack(fill='x', padx=5, pady=2)
        tk.Label(loads_row1, text="Nx:", bg='white').pack(side='left'); tk.Entry(loads_row1, textvariable=self.Nx_var, width=8).pack(side='left', padx=4)
        tk.Label(loads_row1, text="Ny:", bg='white').pack(side='left', padx=(10,0)); tk.Entry(loads_row1, textvariable=self.Ny_var, width=8).pack(side='left', padx=4)
        tk.Label(loads_row1, text="Nxy:", bg='white').pack(side='left', padx=(10,0)); tk.Entry(loads_row1, textvariable=self.Nxy_var, width=8).pack(side='left', padx=4)
        loads_row2 = tk.Frame(self.load_section, bg='white'); loads_row2.pack(fill='x', padx=5, pady=2)
        tk.Label(loads_row2, text="Mx:", bg='white').pack(side='left'); tk.Entry(loads_row2, textvariable=self.Mx_var, width=8).pack(side='left', padx=4)
        tk.Label(loads_row2, text="My:", bg='white').pack(side='left', padx=(10,0)); tk.Entry(loads_row2, textvariable=self.My_var, width=8).pack(side='left', padx=4)
        tk.Label(loads_row2, text="Mxy:", bg='white').pack(side='left', padx=(10,0)); tk.Entry(loads_row2, textvariable=self.Mxy_var, width=8).pack(side='left', padx=4)
        loads_row3 = tk.Frame(self.load_section, bg='white'); loads_row3.pack(fill='x', padx=5, pady=2)
        tk.Label(loads_row3, text="ΔT:", bg='white').pack(side='left'); tk.Entry(loads_row3, textvariable=self.deltaT_var, width=8).pack(side='left', padx=4)
    # (Allowables removed here; now per-lamina in Add New Lamina)

        # --- Optimization options section ---
        opt_chk_frame = tk.Frame(parent, bg='white')
        opt_chk_frame.pack(fill='x', padx=5, pady=(10,0))
        opt_chk = ttk.Checkbutton(opt_chk_frame, text="Optimize Thickness", variable=self.optimize_var, command=self.toggle_opt_section)
        opt_chk.pack(anchor='w')
        self.opt_section = tk.LabelFrame(parent, text="Optimization Options", font=("Arial", 11, "bold"), bg='white')
        # inside opt_section (initially hidden)
        opt_row1 = tk.Frame(self.opt_section, bg='white'); opt_row1.pack(fill='x', padx=5, pady=2)
        tk.Label(opt_row1, text="Desired_SM:", bg='white').pack(side='left'); tk.Entry(opt_row1, textvariable=self.SM_wanted_var, width=8).pack(side='left', padx=4)
        same_chk = ttk.Checkbutton(opt_row1, text="Same thickness per lam", variable=self.same_thickness_var)
        same_chk.pack(side='left', padx=(10,0))
        opt_row2 = tk.Frame(self.opt_section, bg='white'); opt_row2.pack(fill='x', padx=5, pady=2)
        tk.Label(opt_row2, text="Total thickness target (mm):", bg='white').pack(side='left'); tk.Entry(opt_row2, textvariable=self.total_thickness_target_mm_var, width=8).pack(side='left', padx=4)
        opt_row3 = tk.Frame(self.opt_section, bg='white'); opt_row3.pack(fill='x', padx=5, pady=2)
        tk.Label(opt_row3, text="lam min/max t (mm):", bg='white').pack(side='left')
        tk.Entry(opt_row3, textvariable=self.lam_min_thick_mm_var, width=6).pack(side='left', padx=3)
        tk.Entry(opt_row3, textvariable=self.lam_max_thick_mm_var, width=6).pack(side='left', padx=3)
        tk.Label(opt_row3, text="suggested: (0.1, 0.5)", fg='gray', bg='white').pack(side='left', padx=6)
        opt_row4 = tk.Frame(self.opt_section, bg='white'); opt_row4.pack(fill='x', padx=5, pady=2)
        tk.Label(opt_row4, text="Criteria:", bg='white').pack(side='left')
        stress_chk = ttk.Checkbutton(opt_row4, text="Stress", variable=self.stress_criteria_var, command=self.on_stress_toggle)
        stress_chk.pack(side='left', padx=(6,2))
        strain_chk = ttk.Checkbutton(opt_row4, text="Strain", variable=self.strain_criteria_var, command=self.on_strain_toggle)
        strain_chk.pack(side='left')

        # Bottom control buttons: Calculate, Reset, and Export CSV
        bottom_buttons = tk.Frame(parent, bg='white')
        bottom_buttons.pack(fill='x', padx=5, pady=12, side='bottom')
        calc_btn = tk.Button(bottom_buttons, text="CALCULATE", command=self.update_results_display,
                             font=("Arial", 14, "bold"), bg='darkgreen', fg='white',
                             width=22, height=2)
        calc_btn.pack(pady=5)
        reset_btn = tk.Button(bottom_buttons, text="Reset", command=self.reset_values,
                              font=("Arial", 10), bg='orange', fg='white', width=18)
        reset_btn.pack(pady=2)
        export_btn = tk.Button(bottom_buttons, text="Export results as .CSV", command=self.export_results_dialog,
                               font=("Arial", 10), bg='darkblue', fg='white', width=22)
        export_btn.pack(pady=(2,0))
        load_btn = tk.Button(bottom_buttons, text="Load from result", command=self.load_results_dialog,
                               font=("Arial", 10), bg='darkblue', fg='white', width=22)
        load_btn.pack(pady=(2,0))

    def update_laminae_display(self):
        """ Update the display of laminae in the treeview """
        for item in self.laminae_tree.get_children():
            self.laminae_tree.delete(item)
            
        for i, lamina in enumerate(self.laminae_list):
            self.laminae_tree.insert('', 'end', values=(
                i + 1,
                f"{lamina['E1']:.1f}",
                f"{lamina['E2']:.1f}",
                f"{lamina['G12']:.2f}",
                f"{lamina['nu12']:.3f}",
                f"{lamina['thickness']:.3f}",
                f"{lamina['angle']:.1f}",
                f"{lamina.get('alpha_1', self.alpha1_micro_var.get()*1e-6)*1e6:.2f}",
                f"{lamina.get('alpha_2', self.alpha2_micro_var.get()*1e-6)*1e6:.2f}",
                f"{lamina.get('sigma_1_LT', 0.0):.1f}",
                f"{lamina.get('sigma_2_LT', 0.0):.1f}",
                f"{lamina.get('sigma_1_LC', 0.0):.1f}",
                f"{lamina.get('sigma_2_LC', 0.0):.1f}",
                f"{lamina.get('tau_12_L', lamina.get('tao_L', 0.0)):.1f}"
            ))
        self.draw_laminate_visual()
        self.thickness_vector = [lamina['thickness']/1000 for lamina in self.laminae_list] # in meters
            
    def add_lamina(self):
        """ Add a new lamina to the list """
        lamina_data = {
            'E1': self.E1_var.get(),
            'E2': self.E2_var.get(),
            'G12': self.G12_var.get(),
            'nu12': self.nu12_var.get(),
            'thickness': self.thickness_var.get(),
            'angle': self.angle_var.get(),
            # per-lamina allowables in MPa
            'sigma_1_LT': self.sigma_1_LT_MPa.get(),
            'sigma_2_LT': self.sigma_2_LT_MPa.get(),
            'sigma_1_LC': self.sigma_1_LC_MPa.get(),
            'sigma_2_LC': self.sigma_2_LC_MPa.get(),
            'tau_12_L': self.tao_L_MPa.get(),
            # thermal expansion (store real values in 1/°C)
            'alpha_1': self.alpha1_micro_var.get()*1e-6,
            'alpha_2': self.alpha2_micro_var.get()*1e-6,
        }
        self.laminae_list.append(lamina_data)
        self.update_laminae_display()
        
    def remove_lamina(self):
        """ Remove selected lamina """
        selected = self.laminae_tree.selection()
        if selected:
            item = self.laminae_tree.item(selected[0])
            layer_index = int(item['values'][0]) - 1
            if 0 <= layer_index < len(self.laminae_list):
                del self.laminae_list[layer_index]
                self.update_laminae_display()
                
    def clear_laminae(self):
        """ Clear all laminae """
        self.laminae_list.clear()
        self.update_laminae_display()
        
    def move_lamina_up(self):
        """ Move selected lamina up """
        selected = self.laminae_tree.selection()
        if selected:
            item = self.laminae_tree.item(selected[0])
            layer_index = int(item['values'][0]) - 1
            if layer_index > 0:
                self.laminae_list[layer_index], self.laminae_list[layer_index-1] = \
                    self.laminae_list[layer_index-1], self.laminae_list[layer_index]
                self.update_laminae_display()
                
    def move_lamina_down(self):
        """ Move selected lamina down """
        selected = self.laminae_tree.selection()
        if selected:
            item = self.laminae_tree.item(selected[0])
            layer_index = int(item['values'][0]) - 1
            if layer_index < len(self.laminae_list) - 1:
                self.laminae_list[layer_index], self.laminae_list[layer_index+1] = \
                    self.laminae_list[layer_index+1], self.laminae_list[layer_index]
                self.update_laminae_display()

    # --- Clipboard operations for laminae ---
    def _get_selected_indices(self):
        selected = self.laminae_tree.selection()
        idxs = []
        for iid in selected:
            try:
                layer = int(self.laminae_tree.item(iid)['values'][0]) - 1
            except Exception:
                continue
            if 0 <= layer < len(self.laminae_list):
                idxs.append(layer)
        return sorted(set(idxs))

    def copy_laminae(self):
        idxs = self._get_selected_indices()
        if not idxs:
            return
        self._lamina_clipboard = [copy.deepcopy(self.laminae_list[i]) for i in idxs]

    def cut_laminae(self):
        idxs = self._get_selected_indices()
        if not idxs:
            return
        self._lamina_clipboard = [copy.deepcopy(self.laminae_list[i]) for i in idxs]
        for i in sorted(idxs, reverse=True):
            del self.laminae_list[i]
        self.update_laminae_display()

    def paste_laminae(self):
        if not getattr(self, '_lamina_clipboard', None):
            return
        idxs = self._get_selected_indices()
        insert_pos = (max(idxs) + 1) if idxs else len(self.laminae_list)
        for d in self._lamina_clipboard:
            self.laminae_list.insert(insert_pos, copy.deepcopy(d))
            insert_pos += 1
        self.update_laminae_display()
                
    def load_preset_1(self):
        # Load [0/45/-45/90]s preset using current input values
        self.clear_laminae()
        angles = [0, 45, -45, 90, 90, -45, 45, 0]
        for angle in angles:
            lamina_data = {
                'E1': self.E1_var.get(),
                'E2': self.E2_var.get(),
                'G12': self.G12_var.get(),
                'nu12': self.nu12_var.get(),
                'thickness': self.thickness_var.get(),
                'angle': angle,
                # copy current per-lamina allowables into each preset ply (MPa)
                'sigma_1_LT': self.sigma_1_LT_MPa.get(),
                'sigma_2_LT': self.sigma_2_LT_MPa.get(),
                'sigma_1_LC': self.sigma_1_LC_MPa.get(),
                'sigma_2_LC': self.sigma_2_LC_MPa.get(),
                'tau_12_L': self.tao_L_MPa.get(),
                'alpha_1': self.alpha1_micro_var.get()*1e-6,
                'alpha_2': self.alpha2_micro_var.get()*1e-6,
            }
            self.laminae_list.append(lamina_data)
        self.update_laminae_display()
        
    def load_preset_2(self):
        # Load [0/90]s preset using current input values
        self.clear_laminae()
        angles = [0, 90, 90, 0]
        for angle in angles:
            lamina_data = {
                'E1': self.E1_var.get(),
                'E2': self.E2_var.get(),
                'G12': self.G12_var.get(),
                'nu12': self.nu12_var.get(),
                'thickness': self.thickness_var.get(),
                'angle': angle,
                'sigma_1_LT': self.sigma_1_LT_MPa.get(),
                'sigma_2_LT': self.sigma_2_LT_MPa.get(),
                'sigma_1_LC': self.sigma_1_LC_MPa.get(),
                'sigma_2_LC': self.sigma_2_LC_MPa.get(),
                'tau_12_L': self.tao_L_MPa.get(),
                'alpha_1': self.alpha1_micro_var.get()*1e-6,
                'alpha_2': self.alpha2_micro_var.get()*1e-6,
            }
            self.laminae_list.append(lamina_data)
        self.update_laminae_display()
        
    def load_preset_3(self):
        # Load [±45]s preset using current input values
        self.clear_laminae()
        angles = [45, -45, -45, 45]
        for angle in angles:
            lamina_data = {
                'E1': self.E1_var.get(),
                'E2': self.E2_var.get(),
                'G12': self.G12_var.get(),
                'nu12': self.nu12_var.get(),
                'thickness': self.thickness_var.get(),
                'angle': angle,
                'sigma_1_LT': self.sigma_1_LT_MPa.get(),
                'sigma_2_LT': self.sigma_2_LT_MPa.get(),
                'sigma_1_LC': self.sigma_1_LC_MPa.get(),
                'sigma_2_LC': self.sigma_2_LC_MPa.get(),
                'tau_12_L': self.tao_L_MPa.get(),
                'alpha_1': self.alpha1_micro_var.get()*1e-6,
                'alpha_2': self.alpha2_micro_var.get()*1e-6,
            }
            self.laminae_list.append(lamina_data)
        self.update_laminae_display()

    def reset_values(self):
        """Reset new lamina input values (MPa units)."""
        self.E1_var.set(125000.0)
        self.E2_var.set(12500.0)
        self.G12_var.set(6890.0)
        self.nu12_var.set(0.38)
        self.thickness_var.set(0.2)
        self.angle_var.set(0.0)
        self.clear_laminae()

    # --- Simple numeric validation applied after UI creation ---
    def _attach_numeric_validation_all(self):
        def is_valid(s: str) -> bool:
            if s == "":
                return False
            if "," in s:
                return False
            try:
                float(s)
                return True
            except Exception:
                return False
        def bind_entry(e: tk.Entry):
            e._last_valid = e.get()
            def on_focus_out(_evt=None):
                val = e.get().strip()
                if not is_valid(val):
                    messagebox.showerror("Input error", "Use '.' for decimals and only numeric characters (no commas).")
                    e.delete(0, tk.END)
                    e.insert(0, getattr(e, '_last_valid', ''))
                else:
                    e._last_valid = val
            e.bind('<FocusOut>', on_focus_out)
        def recurse(w):
            for child in w.winfo_children():
                if isinstance(child, tk.Entry):
                    bind_entry(child)
                recurse(child)
        recurse(self.input_content)

### Managing the results section and the outputs ###
    def update_results_display(self):
        """ Update the results display in the results tabs """
        self.results = self.calculate_everything(self.thickness_vector)
        Vector_Q_matrices = self.results['Vector_Q_matrices']
        A = self.results['A']
        B = self.results['B']
        D = self.results['D']
        s1 = self.results['s1']
        A2 = self.results['A2']
        A3 = self.results['A3']
        A4 = self.results['A4']
        A5 = self.results['A5']
        Ex = self.results['Ex']
        Ey = self.results['Ey']
        Gxy = self.results['Gxy']
        nuxy = self.results['nuxy']
        nuyx = self.results['nuyx']
        Exb = self.results['Exb']
        Eyb = self.results['Eyb']
        Dx = self.results['Dx']
        Dy = self.results['Dy']
        t = self.results['t']
        self.update_matrices_display(Vector_Q_matrices,A, B, D, s1, A2, A3, A4, A5)
        self.refresh_properties_display(Ex, Ey, Gxy, nuxy, nuyx, Exb, Eyb, Dx, Dy, t, s1)
        if self.add_load_var.get():
            # Assess failure and update related views
            assess = self.assess_failure(self.results)
            # store failed plies map for drawing
            self.failed_plies = {f['index']: f['reasons'] for f in assess['failures']}
            # redraw laminate with failures highlighted
            self.draw_laminate_visual()
            # update checks tab
            self.update_checks_tab(assess)
        # If failed and optimization enabled (and Load Condition active), run optimization
        if self.optimize_var.get() and self.add_load_var.get():
            self.optimize_laminate()
    
    def setup_results_widgets(self):
        """ Setup widgets for the matrices in the matrices tab (matrices_frame) """
        
        matrices_canvas = tk.Canvas(self.matrices_frame, bg='white') #the canvas allows scrolling
        # it is putted inside the tab
        matrices_scrollbar = ttk.Scrollbar(self.matrices_frame, orient="vertical", 
                                         command=matrices_canvas.yview)
        self.matrices_content = tk.Frame(matrices_canvas, bg='white') # this frame holds the actual content
        
        matrices_canvas.configure(yscrollcommand=matrices_scrollbar.set) # configure the canvas to use the scrollbar
        matrices_canvas.pack(side="left", fill="both", expand=True) # pack the canvas to fill the tab
        matrices_scrollbar.pack(side="right", fill="y") # pack the scrollbar to the right side
        matrices_canvas.create_window((0, 0), window=self.matrices_content, anchor="nw") # put the content frame inside the canvas

        # Setup widgets for the properties in the properties tab (properties_frame)
        properties_canvas = tk.Canvas(self.properties_frame, bg='white')
        properties_scrollbar = ttk.Scrollbar(self.properties_frame, orient="vertical", 
                                           command=properties_canvas.yview)
        self.properties_content = tk.Frame(properties_canvas, bg='white')
        
        properties_canvas.configure(yscrollcommand=properties_scrollbar.set)
        properties_canvas.pack(side="left", fill="both", expand=True)
        properties_scrollbar.pack(side="right", fill="y")
        properties_canvas.create_window((0, 0), window=self.properties_content, anchor="nw")

        # Update scroll region
        def configure_scroll_matrices(event):
            matrices_canvas.configure(scrollregion=matrices_canvas.bbox("all"))
        def configure_scroll_properties(event):
            properties_canvas.configure(scrollregion=properties_canvas.bbox("all"))
            
        self.matrices_content.bind("<Configure>", configure_scroll_matrices) # bind the configure event to update the scroll region
        # every time the content frame is resized
        self.properties_content.bind("<Configure>", configure_scroll_properties)

    # Setup for visual representation tab
        self.visual_canvas = tk.Canvas(self.visual_frame, bg='white')
        visual_scrollbar = ttk.Scrollbar(self.visual_frame, orient="vertical", command=self.visual_canvas.yview)
        self.visual_canvas.configure(yscrollcommand=visual_scrollbar.set)
        
        self.visual_canvas.pack(side="left", fill="both", expand=True)
        visual_scrollbar.pack(side="right", fill="y")

        def configure_scroll_visual(event):
            # Pass the correct width from the event to the drawing function
            self.draw_laminate_visual(event.width)
            self.visual_canvas.configure(scrollregion=self.visual_canvas.bbox("all"))
        
        self.visual_canvas.bind("<Configure>", configure_scroll_visual)
        # Setup for checks tab (scrollable)
        self.checks_canvas = tk.Canvas(self.checks_frame, bg='white')
        checks_scrollbar = ttk.Scrollbar(self.checks_frame, orient="vertical", command=self.checks_canvas.yview)
        self.checks_content = tk.Frame(self.checks_canvas, bg='white')
        self.checks_canvas.configure(yscrollcommand=checks_scrollbar.set)
        self.checks_canvas.pack(side="left", fill="both", expand=True)
        checks_scrollbar.pack(side="right", fill="y")
        self.checks_canvas.create_window((0,0), window=self.checks_content, anchor='nw')
        def configure_scroll_checks(event):
            self.checks_canvas.configure(scrollregion=self.checks_canvas.bbox("all"))
        self.checks_content.bind("<Configure>", configure_scroll_checks)
        # Setup for optimization tab (simple frame, will be filled on demand)
        
    def update_matrices_display(self,Vector_Q_matrices, A, B, D, angles=None, A2=None, A3=None, A4=None, A5=None):
        """Method to update the matrices display in the matrices tab"""
        # Clean the frame
        for widget in self.matrices_content.winfo_children():
            widget.destroy()
        # Colors for the matrices
        color_A = '#ffcccc'  # Light red for A matrix
        color_B = '#ccffcc'  # Light green for B matrix  
        color_D = '#ccccff'  # Light blue for D matrix
        # First display the force-displacement relationship
        relationship_title = tk.Label(self.matrices_content, 
                                    text="FORCE-DISPLACEMENT RELATIONSHIP", 
                                    font=("Arial", 16, "bold"), bg='white', fg='navy')
        relationship_title.grid(row=0, column=5, columnspan=1, pady=(10,20))
        # Create the relationship display frame
        rel_frame = tk.Frame(self.matrices_content, bg='white')
        rel_frame.grid(row=1, column=1, columnspan=12, pady=10)
        # Left side - Forces vector (column vector)
        forces_frame = tk.Frame(rel_frame, bg='lightgray', relief='solid', bd=2)
        forces_frame.grid(row=0, column=0, padx=5)
        # Forces vector elements
        force_elements = ["Nx", "Ny", "Nxy", "Mx", "My", "Mxy"]
        for i, element in enumerate(force_elements):
            tk.Label(forces_frame, text=element, font=("Arial", 13, "bold"), bg='lightgray').grid(row=i, column=0, sticky='ew', padx=5, pady=2)
        equal_label = tk.Label(rel_frame, text="=", font=("Arial", 14, "bold"), bg='white')
        equal_label.grid(row=0, column=1, padx=10)
        # Matrix structure with classical [A B; B D] layout
        matrix_frame = tk.Frame(rel_frame, bg='white', relief='solid', bd=2)
        matrix_frame.grid(row=0, column=2, padx=10)
        # A matrix (top-left)
        A_frame = tk.Frame(matrix_frame, bg=color_A, relief='solid', bd=1)
        A_frame.grid(row=0, column=0, padx=2, pady=2)
        for i in range(3):
            for j in range(3):
                value = f"{A[i,j]:.2e}" if abs(A[i,j]) >= 1e-7 else "0"
                tk.Label(A_frame, text=value, font=("Courier", 10), bg=color_A, width=8).grid(row=i, column=j, padx=1, pady=1)
        # B matrix (top-right)
        B_frame = tk.Frame(matrix_frame, bg=color_B, relief='solid', bd=1)
        B_frame.grid(row=0, column=1, padx=2, pady=2)
        for i in range(3):
            for j in range(3):
                value = f"{B[i,j]:.2e}" if abs(B[i,j]) >= 1e-7 else "0"
                tk.Label(B_frame, text=value, font=("Courier", 10), bg=color_B, width=8).grid(row=i, column=j, padx=1, pady=1)
        # B matrix (bottom-left) - same as top-right
        B_frame2 = tk.Frame(matrix_frame, bg=color_B, relief='solid', bd=1)
        B_frame2.grid(row=1, column=0, padx=2, pady=2)
        for i in range(3):
            for j in range(3):
                value = f"{B[i,j]:.2e}" if abs(B[i,j]) >= 1e-7 else "0"
                tk.Label(B_frame2, text=value, font=("Courier", 10), bg=color_B, width=8).grid(row=i, column=j, padx=1, pady=1)
        # D matrix (bottom-right)
        D_frame = tk.Frame(matrix_frame, bg=color_D, relief='solid', bd=1)
        D_frame.grid(row=1, column=1, padx=2, pady=2)
        for i in range(3):
            for j in range(3):
                value = f"{D[i,j]:.2e}" if abs(D[i,j]) >= 1e-7 else "0"
                tk.Label(D_frame, text=value, font=("Courier", 10), bg=color_D, width=8).grid(row=i, column=j, padx=1, pady=1)
        # Right side - Strains vector (column vector)
        strains_frame = tk.Frame(rel_frame, bg='lightyellow', relief='solid', bd=2)
        strains_frame.grid(row=0, column=3, padx=10)
        # Strains vector elements
        strain_elements = ["εx", "εy", "γxy", "κx", "κy", "κxy"]
        for i, element in enumerate(strain_elements):
            tk.Label(strains_frame, text=element, font=("Arial", 13, "bold"), bg='lightyellow').grid(row=i, column=0, sticky='ew', padx=5, pady=2)
        # Legend for classical [A B; B D] layout
        legend_frame = tk.Frame(self.matrices_content, bg='white')
        legend_frame.grid(row=3, column=1, columnspan=12, pady=20)
        tk.Label(legend_frame, text="MATRIX LEGEND :", font=("Arial", 12, "bold"), bg='white').grid(row=0, column=0, columnspan=6, pady=5)
        # A matrix legend
        tk.Label(legend_frame, text="  A  ", font=("Arial", 10, "bold"), bg=color_A, relief='solid', bd=1).grid(row=1, column=0, padx=5)
        tk.Label(legend_frame, text="Extensional Stiffness Matrix (N/m)", font=("Arial", 10), bg='white').grid(row=1, column=1, sticky='w', padx=5)
        # B matrix legend  
        tk.Label(legend_frame, text="  B  ", font=("Arial", 10, "bold"), bg=color_B, relief='solid', bd=1).grid(row=1, column=2, padx=5)
        tk.Label(legend_frame, text="Coupling Stiffness Matrix (N)", font=("Arial", 10), bg='white').grid(row=1, column=3, sticky='w', padx=5)
        # D matrix legend
        tk.Label(legend_frame, text="  D  ", font=("Arial", 10, "bold"), bg=color_D, relief='solid', bd=1).grid(row=1, column=4, padx=5)
        tk.Label(legend_frame, text="Bending Stiffness Matrix (N⋅m)", font=("Arial", 10), bg='white').grid(row=1, column=5, sticky='w', padx=5)
        # Special Case classification
        special_frame = tk.Frame(self.matrices_content, bg='white')
        special_frame.grid(row=2, column=1, columnspan=12, pady=(0,10))
        # Determine special cases if angles are provided
        special_text = "Special Case: "
        if angles is not None:
            special_cases = self.classify_laminate(angles, A, B, D, A2, A3, A4, A5)
            if special_cases:
                special_text += ", ".join(special_cases)
            else:
                special_text += "None"
        else:
            special_text += "None"
        special_label = tk.Label(special_frame, text=special_text, 
                                font=("Arial", 12, "bold"), bg='white', fg='darkgreen')
        special_label.grid(row=0, column=0, pady=5)
        # Separator line before individual matrices
        separator = ttk.Separator(self.matrices_content, orient='horizontal')
        separator.grid(row=4, column=1, columnspan=12, sticky='ew')
        # Individual matrices title
        individual_title = tk.Label(self.matrices_content, 
                                   text="INDIVIDUAL MATRICES", 
                                   font=("Arial", 14, "bold"), bg='white', fg='darkblue')
        individual_title.grid(row=5, column=0, columnspan=12, pady=(10,20))
        # Create a frame for the matrices in 2x2 layout: A-D on top, B centered below
        matrices_display_frame = tk.Frame(self.matrices_content, bg='white')
        matrices_display_frame.grid(row=6, column=0, columnspan=12, pady=10)
        # A matrix (top-left)
        A_container = tk.Frame(matrices_display_frame, bg='white')
        A_container.grid(row=0, column=0, padx=20, pady=10)
        A_title = tk.Label(A_container, text="Matrix A (N/m)", 
                          font=("Arial", 12, "bold"), bg='white', fg='darkblue')
        A_title.grid(row=0, column=0, columnspan=3, pady=(0,10))
        for row in range(3):
            for col in range(3):
                num = A[row,col]
                value = f"{num:.3e}" if abs(num) >= 1e-7 else "0.00"
                cell = tk.Label(A_container, text=value, 
                               font=("Courier", 10), bg=color_A, 
                               relief='solid', bd=1, width=12)
                cell.grid(row=row+1, column=col, padx=1, pady=1)
        # D matrix (top-right)
        D_container = tk.Frame(matrices_display_frame, bg='white')
        D_container.grid(row=0, column=1, padx=20, pady=10)
        D_title = tk.Label(D_container, text="Matrix D (N⋅m)", 
                          font=("Arial", 12, "bold"), bg='white', fg='darkblue')
        D_title.grid(row=0, column=0, columnspan=3, pady=(0,10))
        for row in range(3):
            for col in range(3):
                num = D[row,col]
                value = f"{num:.3e}" if abs(num) >= 1e-7 else "0.00"
                cell = tk.Label(D_container, text=value, 
                               font=("Courier", 10), bg=color_D, 
                               relief='solid', bd=1, width=12)
                cell.grid(row=row+1, column=col, padx=1, pady=1)
        # B matrix (bottom, centered under A and D)
        B_container = tk.Frame(matrices_display_frame, bg='white')
        B_container.grid(row=1, column=0, columnspan=2, pady=(20,10))
        B_title = tk.Label(B_container, text="Matrix B (N)", 
                          font=("Arial", 12, "bold"), bg='white', fg='darkblue')
        B_title.grid(row=0, column=0, columnspan=3, pady=(0,10))
        for row in range(3):
            for col in range(3):
                num = B[row,col]
                value = f"{num:.3e}" if abs(num) >= 1e-7 else "0.00"
                cell = tk.Label(B_container, text=value, 
                               font=("Courier", 10), bg=color_B, 
                               relief='solid', bd=1, width=12)
                cell.grid(row=row+1, column=col, padx=1, pady=1)
        # Show individual Q matrices
        # Add separator and title for individual Q matrices
        separator2 = ttk.Separator(self.matrices_content, orient='horizontal')
        separator2.grid(row=15, column=1, columnspan=12, pady=20, sticky='ew')
        individual_q_title = tk.Label(self.matrices_content, 
                                     text="INDIVIDUAL LAMINA Q MATRICES IN GLOBAL REFERENCE FRAME (Pa)", 
                                     font=("Arial", 14, "bold"), bg='white', fg='darkred')
        individual_q_title.grid(row=16, column=2, columnspan=8)
        # Display Q matrices in a grid (2 per row)
        for i, Q_matrix in enumerate(Vector_Q_matrices):
            row_offset = 17 + (i // 2) * 6
            col_offset = (i % 2) * 6 + 1
            # Matrix title
            lamina_info = self.laminae_list[i]
            title = f"Lamina {i+1}: θ={lamina_info['angle']}°, t={lamina_info['thickness']}mm"
            title_label = tk.Label(self.matrices_content, text=title, 
                                  font=("Arial", 10, "bold"), bg='white', fg='darkred')
            title_label.grid(row=row_offset, column=col_offset, columnspan=3, sticky='w', pady=(10,2))
            # Matrix values
            for row in range(3):
                for col in range(3):
                    num = Q_matrix[row,col]
                    value = f"{num:.2e}"
                    if abs(num) < 1e-7: 
                        value = "0.00"
                    cell = tk.Label(self.matrices_content, text=value, 
                                   font=("Courier", 8), bg='mistyrose', 
                                   relief='solid', width=12)
                    cell.grid(row=row_offset+1+row, column=col_offset+col, padx=1, pady=1, sticky='ew')

    def update_checks_tab(self, assess):
        """Update the Checks & Through-Thickness tab contents"""
        # Clear
        for w in self.checks_content.winfo_children():
            w.destroy()
        # Summary
        title = tk.Label(self.checks_content, text="LAMINATE CHECK SUMMARY", font=("Arial", 16, "bold"), bg='white', fg='navy')
        title.pack(pady=(10,5))
        if assess['passed']:
            ok = tk.Label(self.checks_content, text=f"PASSED ✓  |  Min SM (stress): {assess['min_sm_stress']:.3f}  |  Min SM (strain): {assess['min_sm_strain']:.3f}",
                          font=("Arial", 12, "bold"), bg='white', fg='darkgreen')
            ok.pack(pady=(0,10))
        else:
            fail = tk.Label(self.checks_content, text="FAILED ✗", font=("Arial", 12, "bold"), bg='white', fg='red')
            fail.pack(pady=(0,6))
            # list failures
            list_frame = tk.Frame(self.checks_content, bg='white'); list_frame.pack(fill='x', padx=10, pady=5)
            for f in assess['failures']:
                txt = f"Lamina {f['index']+1}: " + ", ".join([f"{t} {d} (FI={val:.2f})" for (t,d,val) in f['reasons']])
                tk.Label(list_frame, text=txt, bg='white', fg='red', anchor='w', justify='left').pack(fill='x')
        # Through-thickness plot: restore three-panel layout (Laminate | Strain | Stress)
        plot_frame = tk.LabelFrame(self.checks_content, text="Through-thickness: laminate | εx(z) | σx(z)", font=("Arial", 11, "bold"), bg='white')
        plot_frame.pack(fill='x', padx=10, pady=10)
        # Slightly smaller canvas so it fits without horizontal scrolling
        canvas = tk.Canvas(plot_frame, bg='white', height=320, width=820)
        canvas.pack(padx=10, pady=10)

        stresses = np.array(self.results.get('stress_i_global_list', []))  # shape [n,3]
        z_edges = np.array(self.results.get('z', []))  # length nplies+1
        if stresses.size > 0 and z_edges.size >= 2:
            # Geometry
            top = 24; bottom = 280
            H = bottom - top
            # Panel widths (reduced to fit inside tab)
            W_geom = 140
            W_plot = 280
            gap = 20
            # Bounds
            Lg = 50; Rg = Lg + W_geom
            L1 = Rg + gap; R1 = L1 + W_plot  # strain panel
            L2 = R1 + gap; R2 = L2 + W_plot  # stress panel

            # Common y-mapper based on thickness
            zmin, zmax = float(np.min(z_edges)), float(np.max(z_edges))
            def ymap(zv: float) -> float:
                # Inverted orientation so Lamina 1 is on top: zmin (bottom) -> top, zmax (top) -> bottom
                return top + (zv - zmin)/(zmax - zmin + 1e-12)*H

            # Mid-plane line across all panels (no text label)
            y_mid = ymap(0.0)
            canvas.create_line(Lg, y_mid, R2, y_mid, fill='#888', dash=(4,3))

            # Draw ply boundaries as faint lines across all panels
            for ze in z_edges:
                yb = ymap(float(ze))
                canvas.create_line(Lg, yb, R2, yb, fill='#eee')

            # Panel 1: Laminate geometry (stack)
            nplies = len(z_edges) - 1
            # alternating colors for plies
            colors = ['#f0f8ff', '#e6f2ff']
            for i in range(nplies):
                y1 = ymap(float(z_edges[i]))
                y2 = ymap(float(z_edges[i+1]))
                canvas.create_rectangle(Lg, y2, Rg, y1, fill=colors[i % 2], outline='#bbb')
                # angle label centered in ply
                if i < len(self.laminae_list):
                    ang = self.laminae_list[i].get('angle', '')
                    canvas.create_text((Lg+Rg)//2, (y1+y2)/2, text=f"{ang}°", fill='#333', font=("Arial", 9))
            canvas.create_text((Lg+Rg)//2, top-8, text='Laminate', fill='#444', font=("Arial", 10))
            # z-axis with direction arrow in front of laminate (thin cartesian-style)
            y_top_surf = ymap(float(z_edges[-1]))  # physical top surface (now at lower y)
            y_bot_surf = ymap(float(z_edges[0]))   # physical bottom surface (now at upper y)
            xz = Lg - 18
            # draw downward-pointing z-axis: from bottom surface (upper y) down to top surface (lower y)
            canvas.create_line(xz, y_bot_surf, xz, y_top_surf, fill='#555', width=1, arrow=tk.LAST, arrowshape=(7,9,3))

            # Panel 2: Strain εx(z)
            eps0 = np.array(self.results.get('eps0', [0,0,0]), dtype=float)
            kappa = np.array(self.results.get('kappa', [0,0,0]), dtype=float)
            z_samp = np.linspace(zmin, zmax, 120)
            epsx = eps0[0] + z_samp * kappa[0]
            # scale symmetrically
            max_abs_eps = float(np.max(np.abs(epsx)))
            if not np.isfinite(max_abs_eps) or max_abs_eps == 0.0:
                max_abs_eps = 1.0
            def xmap_eps(val: float) -> float:
                return L1 + (val/max_abs_eps)*(W_plot/2) + W_plot/2
            # ε=0 vertical axis with ticks at ply interfaces
            x0_eps = xmap_eps(0.0)
            canvas.create_line(x0_eps, top, x0_eps, bottom, fill='#333')
            for ze in z_edges:
                yb = ymap(float(ze))
                canvas.create_line(x0_eps-4, yb, x0_eps+4, yb, fill='#333')
            # strain polyline with arrows
            pts = [(xmap_eps(float(e)), ymap(float(z))) for e, z in zip(epsx, z_samp)]
            for a, b in zip(pts[:-1], pts[1:]):
                canvas.create_line(a[0], a[1], b[0], b[1], fill='#1f77b4', width=2)
            # arrowheads along the curve
            if len(pts) > 8:
                for idx in np.linspace(6, len(pts)-7, 6, dtype=int):
                    x0, y0 = pts[idx]
                    x1, y1 = pts[idx+1]
                    dx, dy = x1-x0, y1-y0
                    nrm = (dx*dx + dy*dy)**0.5 + 1e-12
                    ux, uy = dx/nrm, dy/nrm
                    size = 9
                    ax, ay = x1 - ux*size, y1 - uy*size
                    px, py = -uy, ux
                    w = 3
                    poly = [x1, y1, ax+px*w, ay+py*w, ax-px*w, ay-py*w]
                    canvas.create_polygon(poly, fill='#1f77b4', outline='')
            # close εx(z) to ε=0 axis at top and bottom with horizontal segments
            e_top = float(eps0[0] + zmax * kappa[0])
            e_bot = float(eps0[0] + zmin * kappa[0])
            xtop = xmap_eps(e_top)
            xbot = xmap_eps(e_bot)
            canvas.create_line(x0_eps, y_top_surf, xtop, y_top_surf, fill='#1f77b4')
            canvas.create_line(x0_eps, y_bot_surf, xbot, y_bot_surf, fill='#1f77b4')
            canvas.create_text((L1+R1)//2, top-8, text='Strain εx(z)', fill='#1f77b4')

            # Panel 3: Stress σx(z)
            # Calculate stress at top and bottom of each ply
            stress_points = [] # List of (sigma_x, z) tuples
            max_sig_val = 0.0
            dT = self.deltaT_var.get()
            
            for i in range(nplies):
                lam = self.laminae_list[i]
                z_bot = z_edges[i]
                z_top = z_edges[i+1]
                
                # Global strains at bot and top
                eps_bot = eps0 + z_bot * kappa
                eps_top = eps0 + z_top * kappa
                
                # Transformation matrix
                angle = lam['angle']
                T_i = self.Theta(angle)
                T_i_transp = T_i.T
                
                # Local strains
                eps_local_bot = T_i_transp @ eps_bot
                eps_local_top = T_i_transp @ eps_top
                
                # Thermal effects
                a1 = lam.get('alpha_1', self.alpha1_micro_var.get()*1e-6)
                a2 = lam.get('alpha_2', self.alpha2_micro_var.get()*1e-6)
                alpha_vec = np.array([a1, a2, 0.0])
                
                eps_local_eff_bot = eps_local_bot - dT * alpha_vec
                eps_local_eff_top = eps_local_top - dT * alpha_vec
                
                # Local Q matrix
                E_1 = lam['E1'] * 1e6
                E_2 = lam['E2'] * 1e6
                G_12 = lam['G12'] * 1e6
                nu_12 = lam['nu12']
                Q_11 = E_1/(1-(nu_12**2)*(E_2/E_1))
                Q_12 = (nu_12*E_2)/(1-(nu_12**2)*(E_2/E_1))
                Q_22 = E_2/(1-(nu_12**2)*(E_2/E_1))
                Q_66 = G_12
                Q_local = np.array([[Q_11, Q_12, 0],[Q_12, Q_22, 0],[0, 0, Q_66]])
                
                # Local stresses
                sig_local_bot = Q_local @ eps_local_eff_bot
                sig_local_top = Q_local @ eps_local_eff_top
                
                # Global stresses
                sig_global_bot = T_i @ sig_local_bot
                sig_global_top = T_i @ sig_local_top
                
                # We are interested in sigma_x (index 0)
                sx_bot = sig_global_bot[0]
                sx_top = sig_global_top[0]
                
                stress_points.append({'z_bot': z_bot, 'z_top': z_top, 'sx_bot': sx_bot, 'sx_top': sx_top})
                max_sig_val = max(max_sig_val, abs(sx_bot), abs(sx_top))

            if max_sig_val == 0.0: max_sig_val = 1.0
            
            def xmap_sig(val: float) -> float:
                return L2 + (val/max_sig_val)*(W_plot/2) + W_plot/2
                
            # σ=0 vertical axis with ticks at ply interfaces
            x0_sig = xmap_sig(0.0)
            canvas.create_line(x0_sig, top, x0_sig, bottom, fill='#333')
            for ze in z_edges:
                yb = ymap(float(ze))
                canvas.create_line(x0_sig-4, yb, x0_sig+4, yb, fill='#333')
                
            # Draw stress profile
            # Close to zero at top
            top_ply_idx = nplies - 1
            if top_ply_idx >= 0:
                top_ply = stress_points[top_ply_idx]
                y_top_global = ymap(top_ply['z_top'])
                x_top_global = xmap_sig(top_ply['sx_top'])
                
                # Line from axis to first point
                canvas.create_line(x0_sig, y_top_global, x_top_global, y_top_global, fill='#d62728', width=3)
                
                for i in range(nplies - 1, -1, -1):
                    ply = stress_points[i]
                    y_t = ymap(ply['z_top'])
                    y_b = ymap(ply['z_bot'])
                    x_t = xmap_sig(ply['sx_top'])
                    x_b = xmap_sig(ply['sx_bot'])
                    
                    # Draw stress variation within ply
                    canvas.create_line(x_t, y_t, x_b, y_b, fill='#d62728', width=3)
                    
                    # If not the last ply (bottom), draw horizontal connection to next ply
                    if i > 0:
                        next_ply = stress_points[i-1] # ply below
                        x_next_top = xmap_sig(next_ply['sx_top'])
                        # y_b is same as y_next_top (z_edges[i])
                        canvas.create_line(x_b, y_b, x_next_top, y_b, fill='#d62728', width=3)
                
                # Close to zero at bottom
                bot_ply = stress_points[0]
                y_bot_global = ymap(bot_ply['z_bot'])
                x_bot_global = xmap_sig(bot_ply['sx_bot'])
                canvas.create_line(x_bot_global, y_bot_global, x0_sig, y_bot_global, fill='#d62728', width=3)

            canvas.create_text((L2+R2)//2, top-8, text='Stress σx(z)', fill='#d62728')

            # removed extra 'z' labels between panels for a cleaner look
        # eps0 and kappa
        ek_frame = tk.Frame(self.checks_content, bg='white'); ek_frame.pack(fill='x', padx=10, pady=10)
        def small_column(frame, vec, title):
            cont = tk.LabelFrame(frame, text=title, font=("Arial", 11, "bold"), bg='white')
            cont.pack(side='left', padx=10)
            for i,v in enumerate(vec):
                tk.Label(cont, text=f"{v:.3e}", font=("Courier", 10), bg='white', width=16, relief='solid', bd=1).grid(row=i, column=0, padx=2, pady=2)
        eps0=self.results['eps0']
        kappa=self.results['kappa']
        for i in range(3):
            if abs(eps0[i])<1e-12: eps0[i]=0.0
            if abs(kappa[i])<1e-12: kappa[i]=0.0
        small_column(ek_frame, eps0, "ε0")
        small_column(ek_frame, kappa, "κ")
        # epsilon and sigma per lamina (global)
        gl_frame = tk.LabelFrame(self.checks_content, text="Per-ply ε and σ (Pa) (global)", font=("Arial", 11, "bold"), bg='white')
        gl_frame.pack(fill='x', padx=10, pady=10)
        for i,(eps, sig) in enumerate(zip(self.results['eps_i_global_list'], self.results['stress_i_global_list'])):
            row = tk.Frame(gl_frame, bg='white'); row.pack(fill='x', padx=5, pady=3)
            tk.Label(row, text=f"Lamina {i+1}", bg='white', width=10).pack(side='left')
            for vec, name, col in [(eps,'ε','#e8f5e9'),(sig,'σ','#fdecea')]:
                box = tk.Frame(row, bg='white'); box.pack(side='left', padx=10)
                tk.Label(box, text=name, bg='white').grid(row=0, column=0, columnspan=1)
                for j,v in enumerate(vec):
                    tk.Label(box, text=f"{v:.3e}", font=("Courier", 9), bg=col, relief='solid', bd=1, width=14).grid(row=j+1, column=0, padx=1, pady=1)

    def update_optimization_tab(self):
        # Clear
        for w in self.optim_frame.winfo_children():
            w.destroy()
        # If no optimization data, show only one informative line
        if not self.optimization_result:
            tk.Label(self.optim_frame, text="No optimization data available", bg='white', font=("Arial", 12)).pack(pady=20)
            return
        title = tk.Label(self.optim_frame, text="OPTIMIZATION RESULTS", font=("Arial", 16, "bold"), bg='white', fg='navy')
        title.pack(pady=(10,10))
        R = self.optimization_result
        info = tk.Frame(self.optim_frame, bg='white'); info.pack(fill='x', padx=10, pady=5)
        x_m = np.array(R['x']) if 'x' in R else np.array([])
        x_mm = x_m*1000.0 if x_m.size>0 else x_m
        total_after_m = float(R.get('total_thickness_after', float('nan')))
        total_after_mm = total_after_m*1000.0 if np.isfinite(total_after_m) else float('nan')
        lines = [
            f"Converged: {R['success']}",
            f"Message: {R['message']}",
            f"Iterations: {R['nit']}",
            f"Objective: {R['fun']:.6f}",
            f"Min SM after optimization: {R.get('min_SM_after', float('nan')):.3f}",
            f"Total thickness after optimization (mm): {total_after_mm:.3f}",
        ]
        for tline in lines:
            tk.Label(info, text=tline, bg='white', anchor='w').pack(fill='x')
        # thickness vector in mm only
        tk.Label(self.optim_frame, text="Optimized thicknesses (mm):", bg='white', font=("Arial", 11, "bold")).pack(anchor='w', padx=10, pady=(10,2))
        tv_mm = tk.Text(self.optim_frame, height=3, bg='#f7f7f7')
        tv_mm.pack(fill='x', padx=10)
        if x_mm.size>0:
            tv_mm.insert('1.0', np.array2string(x_mm, precision=3, separator=', '))
        tv_mm.config(state='disabled')
        # Removed meters view to keep UI consistent in mm
        
        # Conditional note when at max bounds but still negative SM
        try:
            max_bound_m = self.lam_max_thick_mm_var.get()/1000.0
            at_max = (x_m.size>0) and np.allclose(x_m, max_bound_m, rtol=0, atol=1e-9)
        except Exception:
            at_max = False
        if bool(R.get('success', False)) and np.isfinite(total_after_m):
            min_sm_after = float(R.get('min_SM_after', float('nan')))
            if np.isfinite(min_sm_after) and (min_sm_after < 0) and at_max:
                note = (
                    "With these bounds and constraints this is the best achievable Safety Margin "
                    "for this laminate and this load, even if negative."
                )
                tk.Label(self.optim_frame, text=note, bg='white', fg='red', wraplength=900, justify='left').pack(fill='x', padx=10, pady=(8,2))

        # Button to apply optimized thicknesses to current laminate
        btn_frame = tk.Frame(self.optim_frame, bg='white')
        btn_frame.pack(fill='x', padx=10, pady=(12,6))
        tk.Button(btn_frame, text="Use optimized thicknesses as current", command=self.apply_optimized_thicknesses,
                  bg='green', fg='white').pack(side='left')
      
    def refresh_properties_display(self, Ex, Ey, Gxy, nuxy, nuyx, Exb, Eyb, Dx, Dy, t, s1):
        """Method to update the properties display in the properties tab"""
        # Clean the frame
        for widget in self.properties_content.winfo_children():
            widget.destroy()
        # Title
        title_label = tk.Label(self.properties_content, 
                              text="Laminate Engineering Properties", 
                              font=("Arial", 20, "bold"), bg='white', fg='navy')
        title_label.pack(pady=10)
        # Sequence display
        seq_label = tk.Label(self.properties_content, text=f"Sequence: {s1}", 
                            font=("Arial", 16), bg='white', fg='darkblue')
        seq_label.pack(pady=5)
        # Properties table
        tree_frame = tk.Frame(self.properties_content, bg='white')
        tree_frame.pack(pady=20, padx=20, fill='both', expand=True)
        columns = ('property', 'value', 'unit')
        tree = ttk.Treeview(tree_frame, columns=columns, show='headings', height=16)
        # a treeview is like a table to display data in rows and columns
        # Configure style for larger font
        style = ttk.Style()
        style.configure("Properties.Treeview", font=("Arial", 14))
        style.configure("Properties.Treeview.Heading", font=("Arial", 16, "bold"))
        tree.configure(style="Properties.Treeview")
        tree.heading('property', text='Property')
        tree.heading('value', text='Value')
        tree.heading('unit', text='Unit')
        tree.column('property', width=350, anchor='w')
        tree.column('value', width=200, anchor='center')
        tree.column('unit', width=150, anchor='center')
        property_data = [
            ("", "", ""),
            ("Total Thickness", f"{t:.3e}", "m"),
            ("", "", ""),
            ("Ex Module", f"{Ex:.3e}", "Pa"),
            ("Ey Module", f"{Ey:.3e}", "Pa"),
            ("Gxy Module", f"{Gxy:.3e}", "Pa"),
            ("", "", ""),
            ("Poisson vxy", f"{nuxy:.3e}", "-"),
            ("Poisson vyx", f"{nuyx:.3e}", "-"),
            ("", "", ""),
            ("Exb Module", f"{Exb:.3e}", "Pa"),
            ("Eyb Module", f"{Eyb:.3e}", "Pa"),
            ("", "", ""),
            ("Stiffness Dx", f"{Dx:.3e}", "N⋅m"),
            ("Stiffness Dy", f"{Dy:.3e}", "N⋅m"),
        ]
        for prop, val, unit in property_data:
            tree.insert('', 'end', values=(prop, val, unit))
        tree.pack(fill='both', expand=True)

### Visual representation of the laminate stack ###
    def draw_laminate_visual(self, canvas_width=None):
        """Draw the laminate stack on the visual canvas"""
        # Clear canvas
        self.visual_canvas.delete("all")
        # Check if there are laminae to draw
        if not self.laminae_list:
            self.visual_canvas.create_text(300, 200, text="No laminae in the stack.",
                                           font=("Arial", 14), fill='gray')
            return
        # Draw global reference frame in top-left corner
        ref_x, ref_y = 30, 50  # Position for the reference system
        axis_len = 25
        # X-axis
        self.visual_canvas.create_line(ref_x, ref_y, ref_x + axis_len, ref_y, 
                                     arrow=tk.LAST, fill='red', width=2)
        self.visual_canvas.create_text(ref_x + axis_len + 5, ref_y, text='x', 
                                     anchor='w', fill='red', font=("Arial", 10, "bold"))
        # Y-axis
        self.visual_canvas.create_line(ref_x, ref_y, ref_x, ref_y - axis_len, 
                                     arrow=tk.LAST, fill='red', width=2)
        self.visual_canvas.create_text(ref_x, ref_y - axis_len - 5, text='y', 
                                     anchor='s', fill='red', font=("Arial", 10, "bold"))
        # Label for the reference frame
        self.visual_canvas.create_text(ref_x, ref_y + 15, text="Global Frame", 
                                     anchor='w', fill='red', font=("Arial", 10, "italic"))
        # Drawing parameters
        if canvas_width is None or canvas_width < 2:
            canvas_width = self.visual_canvas.winfo_width()
            if canvas_width < 2: return # Don't draw if canvas is not ready
        ply_width = 380 
        ply_height = 225
        skew_offset = 50 # Horizontal skew for the parallelogram
        vertical_gap = 80 # Vertical gap between plies
        
        total_width_of_ply = ply_width + skew_offset
        total_drawing_width = total_width_of_ply 
        start_x = (canvas_width - total_drawing_width) / 2
        if start_x < 10: start_x = 10 # Add a small margin
        start_y = 20 # Moved the whole drawing up
        # Draw from top to bottom
        for i, lamina in enumerate(self.laminae_list):
            y_pos = start_y + i * vertical_gap
            # Define parallelogram points
            p1 = (start_x, y_pos + ply_height)
            p2 = (start_x + ply_width, y_pos + ply_height)
            p3 = (start_x + ply_width + skew_offset, y_pos)
            p4 = (start_x + skew_offset, y_pos)
            angle = lamina['angle']
            color = self.get_color_for_angle(angle)
            failed = i in self.failed_plies
            if failed:
                color = '#ff8a80'
            # Draw ply outline
            self.visual_canvas.create_polygon(p1, p2, p3, p4, fill=color, outline='black', width=1)
            # Draw fiber orientation lines
            # --- Special handling for 0 and 180 degrees ---
            if angle % 180 == 0:
                num_lines = 15
                line_spacing = ply_height / num_lines
                for j in range(1, num_lines):
                    # y position for the horizontal line
                    line_y = y_pos + j * line_spacing
                    # Calculate x coordinates on the skewed edges for this y
                    x_left = p4[0] + (p1[0] - p4[0]) * (line_y - p4[1]) / (p1[1] - p4[1])
                    x_right = p3[0] + (p2[0] - p3[0]) * (line_y - p3[1]) / (p2[1] - p3[1])
                    self.visual_canvas.create_line(x_left, line_y, x_right, line_y, fill='black', width=1)
            else:
                # --- General case for other angles ---
                angle_rad = np.deg2rad(angle)
                center_x = start_x + ply_width / 2 + skew_offset / 2
                center_y = y_pos + ply_height / 2
                # Make line_length large enough to cover the ply at any angle
                line_length = (ply_width + skew_offset + ply_height) * 1.5
                num_lines = 30
                line_spacing = 10 # Increased spacing for better visibility
                cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
                # Define the clipping polygon for this lamina
                clip_polygon = [p1, p2, p3, p4]
                for j in range(-num_lines, num_lines + 1):
                    # Create long lines that will be clipped
                    x_start_local = -line_length / 2
                    x_end_local = line_length / 2
                    y_local = j * line_spacing
                    # Rotate points
                    x1_rot = x_start_local * cos_a - y_local * sin_a
                    y1_rot = x_start_local * sin_a + y_local * cos_a
                    x2_rot = x_end_local * cos_a - y_local * sin_a
                    y2_rot = x_end_local * sin_a + y_local * cos_a
                    # Translate to ply position
                    x1 = center_x + x1_rot
                    y1 = center_y - y1_rot
                    x2 = center_x + x2_rot
                    y2 = center_y - y2_rot
                    # Clip the line against the polygon
                    clipped_line = self.clip_line_to_polygon((x1, y1, x2, y2), clip_polygon)
                    if clipped_line:
                        self.visual_canvas.create_line(clipped_line, fill='black', width=1)
            # Redraw the polygon outline to ensure it's on top
            self.visual_canvas.create_polygon(p1, p2, p3, p4, fill='', outline='black', width=1.5)
            # Draw angle text
            text_x = p3[0] + 10
            text_y = p3[1] + (p1[1] - p4[1]) / 2 - 80 # Moved labels up
            self.visual_canvas.create_text(text_x, text_y, text=f"{lamina['angle']}°",
                                           font=("Arial", 12, "bold"), anchor='w', fill=color)
            # Failure badge
            if failed:
                reasons = self.failed_plies[i]
                tag = "FAIL: " + ", ".join(sorted(set([('σ' if t=='stress' else 'ε')+d for t,d,_ in reasons])))
                self.visual_canvas.create_text(text_x, text_y+18, text=tag, font=("Arial", 10, "bold"), anchor='w', fill='red')
        # Update scroll region
        self.visual_canvas.config(scrollregion=self.visual_canvas.bbox("all"))

    def get_color_for_angle(self, angle):
        """Generates a bright color based on the angle."""
        # Normalize angle to 0-1 range (e.g., for hue)
        # Using a range of -90 to 90 for angles
        normalized_angle = (angle + 90) / 180.0
        # Use HSV color space for bright colors: hue from angle, saturation and value are high
        # Hue is wrapped around to avoid similar red/pinks at both ends
        hue = normalized_angle * 0.8  # Use 80% of the color wheel
        saturation = 0.9
        value = 0.95
        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
        # Convert RGB from 0-1 to 0-255 and format as hex string
        return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'

    def clip_line_to_polygon(self, line, polygon):
        """Clips a line to a convex polygon using Cyrus-Beck algorithm."""
        x1, y1, x2, y2 = line
        t0, t1 = 0.0, 1.0
        line_dx, line_dy = x2 - x1, y2 - y1
        for i in range(len(polygon)):
            p1 = polygon[i]
            p2 = polygon[(i + 1) % len(polygon)]
            # Edge normal (pointing inwards for a clockwise polygon)
            # Our polygon is defined p1,p2,p3,p4 which is clockwise
            edge_dx, edge_dy = p2[0] - p1[0], p2[1] - p1[1]
            normal = (edge_dy, -edge_dx)
            # Vector from edge point to line start
            w = (x1 - p1[0], y1 - p1[1])
            numerator = -(normal[0] * w[0] + normal[1] * w[1]) # -n . w
            denominator = normal[0] * line_dx + normal[1] * line_dy # n . (line_end - line_start)
            if abs(denominator) < 1e-9: # Line is parallel to edge
                if numerator < 0: # Line is outside the edge, so it's clipped
                    return None
                # else, line is inside or on the edge, continue
            else:
                t = numerator / denominator
                if denominator > 0: # Line is entering across this edge
                    t0 = max(t0, t)
                else: # Line is exiting across this edge
                    t1 = min(t1, t)
        if t0 > t1:
            return None # Line is completely outside
        # Calculate clipped line points
        clipped_x1 = x1 + t0 * line_dx
        clipped_y1 = y1 + t0 * line_dy
        clipped_x2 = x1 + t1 * line_dx
        clipped_y2 = y1 + t1 * line_dy
        return clipped_x1, clipped_y1, clipped_x2, clipped_y2

    # --- Helpers for new inputs ---
    def toggle_load_section(self):
        if self.add_load_var.get():
            # Turning ON Loads: show loads section
            self.load_section.pack(fill='x', padx=5, pady=5, side='top')
        else:
            # Turning OFF Loads: hide loads and also force Optimization OFF
            self.load_section.pack_forget()
            if self.optimize_var.get():
                self.optimize_var.set(False)
                # Hide optimization section if visible
                try:
                    self.opt_section.pack_forget()
                except Exception:
                    pass

    def toggle_opt_section(self):
        if self.optimize_var.get():
            # Auto-enable Load Condition when Optimization is enabled
            if not self.add_load_var.get():
                self.add_load_var.set(True)
            # Ensure load section is visible
            self.toggle_load_section()
            # Reorder sections so LOAD always appears BEFORE OPTIMIZE
            try:
                self.opt_section.pack_forget()
            except Exception:
                pass
            try:
                self.load_section.pack_forget()
            except Exception:
                pass
            # Pack Load first, then Optimize underneath
            self.load_section.pack(fill='x', padx=5, pady=5, side='top')
            self.opt_section.pack(fill='x', padx=5, pady=5, side='top')
        else:
            self.opt_section.pack_forget()

    def on_stress_toggle(self):
        """Ensure Stress/Strain checkboxes are mutually exclusive and one is always on."""
        if self.stress_criteria_var.get():
            # Stress turned on -> turn Strain off
            self.strain_criteria_var.set(False)
        else:
            # Prevent both off -> turn Strain on
            self.strain_criteria_var.set(True)

    def on_strain_toggle(self):
        """Ensure Stress/Strain checkboxes are mutually exclusive and one is always on."""
        if self.strain_criteria_var.get():
            # Strain turned on -> turn Stress off
            self.stress_criteria_var.set(False)
        else:
            # Prevent both off -> turn Stress on
            self.stress_criteria_var.set(True)

    def get_loads_and_allowables(self):
        """Read current GUI values and return (N_vec, M_vec, allowables placeholder).

        Note: allowables are now per-lamina and stored within each lamina entry (in MPa),
        so this function returns an empty dict for allowables for backward compatibility.
        """
        N_vec = np.array([self.Nx_var.get(), self.Ny_var.get(), self.Nxy_var.get()], dtype=float)
        M_vec = np.array([self.Mx_var.get(), self.My_var.get(), self.Mxy_var.get()], dtype=float)
        allow = {}
        return N_vec, M_vec, allow

    # --- Credits tab setup and Export helpers ---
    def _setup_credits_tab(self):
        """Populate the Credits tab with copyable text."""
        wrapper = tk.Frame(self.credits_frame, bg='white')
        wrapper.pack(fill='both', expand=True, padx=12, pady=12)
        tk.Label(wrapper, text="Credits", bg='white', fg='navy', font=("Arial", 18, "bold")).pack(anchor='w', pady=(0,8))
        txt = tk.Text(wrapper, height=10, wrap='word', bg='#f7f7f7', font=("Arial", 13))
        txt.pack(fill='both', expand=True)
        content = (
            "Credits of the app go to Francesco Sessa\n"
            "Linkedin Profile: www.linkedin.com/in/francesco-sessa-aer\n"
            "Youtube Profile: www.youtube.com/@francescosessa8861\n"
            "Github Profile: www.github.com/FranciSessa\n"
        )
        txt.insert('1.0', content)
        # Make text copyable but not editable
        txt.config(state='disabled')

    def gather_export_data(self):
        """Collect inputs, outputs, and metadata needed for export."""
        data = {}
        # Credits (fixed)
        data['credits'] = {
            'line1': 'Credits of the app go to Francesco Sessa',
            'linkedin': 'www.linkedin.com/in/francesco-sessa-aer',
            'youtube': 'www.youtube.com/@francescosessa8861',
            'github': 'www.github.com/FranciSessa',
        }
        # Inputs
        data['inputs'] = {
            'laminae': [dict(l) for l in self.laminae_list],
            'thickness_vector_m': list(self.thickness_vector),
            'thickness_vector_mm': [t*1000.0 for t in self.thickness_vector],
            'loads_enabled': bool(self.add_load_var.get()),
            'loads_N_per_m': {
                'Nx': self.Nx_var.get(), 'Ny': self.Ny_var.get(), 'Nxy': self.Nxy_var.get()
            },
            'moments_N_per_m': {
                'Mx': self.Mx_var.get(), 'My': self.My_var.get(), 'Mxy': self.Mxy_var.get()
            },
            'deltaT_degC': self.deltaT_var.get(),
            # per-lamina allowables are embedded within 'laminae'
            'optimization': {
                'enabled': bool(self.optimize_var.get()),
                'SM_wanted': self.SM_wanted_var.get(),
                'same_thickness_per_lam': bool(self.same_thickness_var.get()),
                'total_thickness_target_mm': self.total_thickness_target_mm_var.get(),
                'lam_min_thick_mm': self.lam_min_thick_mm_var.get(),
                'lam_max_thick_mm': self.lam_max_thick_mm_var.get(),
                'criteria': 'stress' if self.stress_criteria_var.get() else 'strain',
            }
        }
        # Results (convert numpy to lists)
        res = self.results if getattr(self, 'results', None) else None
        if res is not None:
            def arr(x):
                return np.array(x).tolist()
            safe_res = {}
            for k,v in res.items():
                try:
                    safe_res[k] = arr(v)
                except Exception:
                    safe_res[k] = v
            data['results'] = safe_res
            # Assess summary
            try:
                data['assess'] = self.assess_failure(res)
            except Exception:
                data['assess'] = None
        else:
            data['results'] = None
            data['assess'] = None
        # Optimization result
        opt = getattr(self, 'optimization_result', None)
        if opt is not None:
            opt_copy = dict(opt)
            try:
                x_m = np.array(opt_copy.get('x', []))
                opt_copy['x_mm'] = (x_m*1000.0).tolist()
            except Exception:
                opt_copy['x_mm'] = []
            data['optimization_result'] = opt_copy
        else:
            data['optimization_result'] = None
        return data

    def load_results_dialog(self):
        """Open a dialog to load inputs from a CSV result file."""
        filename = filedialog.askopenfilename(
            title="Select Result CSV File",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not filename:
            return
        
        try:
            self.import_from_csv(filename)
            tk.messagebox.showinfo("Load Success", "Data loaded successfully from CSV.")
            # Refresh UI
            self.update_laminae_display()
            # Optionally update results if we want to show them immediately
            # self.update_results_display() 
        except Exception as e:
            tk.messagebox.showerror("Load Error", f"Failed to load file:\n{str(e)}")

    def import_from_csv(self, filename):
        import csv
        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            rows = list(reader)
        
        # Helper to find section index
        def find_row_index(start_str):
            for i, r in enumerate(rows):
                if r and start_str in r[0]:
                    return i
            return -1

        # 1. Parse Laminae
        # Look for header: "Laminae (index..."
        lam_header_idx = find_row_index("Laminae (index")
        if lam_header_idx != -1:
            col_header_idx = lam_header_idx + 1
            if col_header_idx < len(rows):
                headers = rows[col_header_idx]
                
                new_laminae = []
                i = col_header_idx + 1
                while i < len(rows):
                    row = rows[i]
                    if not row or not row[0].strip(): # Empty line ends section
                        break
                    
                    lam_data = {}
                    # Helper to get float
                    def get_val(col_name, default=0.0):
                        try:
                            idx = headers.index(col_name)
                            return float(row[idx])
                        except (ValueError, IndexError):
                            return default

                    lam_data['E1'] = get_val("E1_MPa")
                    lam_data['E2'] = get_val("E2_MPa")
                    lam_data['G12'] = get_val("G12_MPa")
                    lam_data['nu12'] = get_val("nu12")
                    lam_data['thickness'] = get_val("thickness_mm")
                    lam_data['angle'] = get_val("angle_deg")
                    
                    # Thermal
                    lam_data['alpha_1'] = get_val("alpha_1_per_degC")
                    lam_data['alpha_2'] = get_val("alpha_2_per_degC")
                    
                    # Allowables
                    lam_data['sigma_1_LT'] = get_val("sigma_1_LT_MPa")
                    lam_data['sigma_2_LT'] = get_val("sigma_2_LT_MPa")
                    lam_data['sigma_1_LC'] = get_val("sigma_1_LC_MPa")
                    lam_data['sigma_2_LC'] = get_val("sigma_2_LC_MPa")
                    
                    # tau_12_L might be named "tau_12_L_MPa"
                    lam_data['tau_12_L'] = get_val("tau_12_L_MPa")
                    
                    new_laminae.append(lam_data)
                    i += 1
                
                if new_laminae:
                    self.laminae_list = new_laminae

        # 2. Parse Loads
        loads_idx = find_row_index("Loads and Moments (N/m)")
        if loads_idx != -1:
            i = loads_idx + 1
            while i < len(rows):
                row = rows[i]
                if not row or not row[0].strip():
                    break
                key = row[0]
                val_str = row[1] if len(row) > 1 else ""
                
                if key == "loads_enabled":
                    self.add_load_var.set(val_str == "True")
                elif key == "Nx": self.Nx_var.set(float(val_str) if val_str else 0.0)
                elif key == "Ny": self.Ny_var.set(float(val_str) if val_str else 0.0)
                elif key == "Nxy": self.Nxy_var.set(float(val_str) if val_str else 0.0)
                elif key == "Mx": self.Mx_var.set(float(val_str) if val_str else 0.0)
                elif key == "My": self.My_var.set(float(val_str) if val_str else 0.0)
                elif key == "Mxy": self.Mxy_var.set(float(val_str) if val_str else 0.0)
                elif key == "DeltaT_degC": self.deltaT_var.set(float(val_str) if val_str else 0.0)
                
                i += 1
            # Update visibility of load section
            self.toggle_load_section()

        # 3. Parse Optimization Options
        opt_idx = find_row_index("Optimization options")
        if opt_idx != -1:
            i = opt_idx + 1
            while i < len(rows):
                row = rows[i]
                if not row or not row[0].strip():
                    break
                key = row[0]
                val_str = row[1] if len(row) > 1 else ""
                
                if key == "enabled":
                    self.optimize_var.set(val_str == "True")
                elif key == "SM_wanted": self.SM_wanted_var.set(float(val_str) if val_str else 0.0)
                elif key == "same_thickness_per_lam": self.same_thickness_var.set(val_str == "True")
                elif key == "total_thickness_target_mm": self.total_thickness_target_mm_var.set(float(val_str) if val_str else 0.0)
                elif key == "lam_min_thick_mm": self.lam_min_thick_mm_var.set(float(val_str) if val_str else 0.0)
                elif key == "lam_max_thick_mm": self.lam_max_thick_mm_var.set(float(val_str) if val_str else 0.0)
                elif key == "criteria":
                    if val_str == "stress":
                        self.stress_criteria_var.set(True)
                        self.strain_criteria_var.set(False)
                    elif val_str == "strain":
                        self.stress_criteria_var.set(False)
                        self.strain_criteria_var.set(True)
                
                i += 1
            # Update visibility of opt section
            self.toggle_opt_section()
    def export_iclc_results(self, path, data):
        """Write a comprehensive CSV export for ICLC results.

        CSV layout (sectioned):
        - Title and credits at the very top
        - Inputs (laminae, thickness vectors, loads, allowables, optimization)
        - Results (engineering constants, ABD, z, eps0, kappa, per-ply lists, FIs)
        - Assessment summary
        - Optimization outcome (if present)
        """
        credits = data.get('credits', {})
        inputs = data.get('inputs', {})
        results = data.get('results')
        assess = data.get('assess')
        optim = data.get('optimization_result')

        def w(writer, row):
            writer.writerow(row)

        def w_blank(writer, n=1):
            for _ in range(n):
                writer.writerow([])

        with open(path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)

            # Header title and credits
            w(writer, ["ICLC output"])  # requested title at the very top
            w(writer, [f"Generated at: {datetime.now().isoformat(timespec='seconds')}"])
            w(writer, [])
            w(writer, ["Credits"])  # simple header
            if credits:
                for k in [
                    'line1', 'linkedin', 'youtube', 'github'
                ]:
                    if k in credits and credits[k]:
                        w(writer, [k, credits[k]])
            else:
                w(writer, ["(no credits provided)"])

            w_blank(writer)
            # Inputs section
            w(writer, ["Inputs"])  # Section header
            # Laminae table (now includes per-lamina allowables in MPa)
            laminae = inputs.get('laminae', [])
            w(writer, ["Laminae (index, E1[MPa], E2[MPa], G12[MPa], nu12, thickness[mm], angle[deg], alpha_1[1/°C], alpha_2[1/°C], σ1_LT[MPa], σ2_LT[MPa], σ1_LC[MPa], σ2_LC[MPa], τ12_L[MPa])"])
            w(writer, ["index", "E1_MPa", "E2_MPa", "G12_MPa", "nu12", "thickness_mm", "angle_deg", "alpha_1_per_degC", "alpha_2_per_degC", "sigma_1_LT_MPa", "sigma_2_LT_MPa", "sigma_1_LC_MPa", "sigma_2_LC_MPa", "tau_12_L_MPa"])
            for i, lam in enumerate(laminae, start=1):
                w(writer, [
                    i,
                    lam.get('E1', ''),
                    lam.get('E2', ''),
                    lam.get('G12', ''),
                    lam.get('nu12', ''),
                    lam.get('thickness', ''),
                    lam.get('angle', ''),
                    lam.get('alpha_1', ''),
                    lam.get('alpha_2', ''),
                    lam.get('sigma_1_LT', ''),
                    lam.get('sigma_2_LT', ''),
                    lam.get('sigma_1_LC', ''),
                    lam.get('sigma_2_LC', ''),
                    lam.get('tau_12_L', lam.get('tao_L', '')),
                ])

            w_blank(writer)
            # Thickness vectors
            w(writer, ["Thickness vectors"])
            w(writer, ["thickness_vector_m"] + list(map(str, inputs.get('thickness_vector_m', []))))
            w(writer, ["thickness_vector_mm"] + list(map(str, inputs.get('thickness_vector_mm', []))))

            w_blank(writer)
            # Loads and moments
            w(writer, ["Loads and Moments (N/m)"])
            loads = inputs.get('loads_N_per_m', {})
            moments = inputs.get('moments_N_per_m', {})
            w(writer, ["loads_enabled", inputs.get('loads_enabled', False)])
            w(writer, ["Nx", loads.get('Nx', '')])
            w(writer, ["Ny", loads.get('Ny', '')])
            w(writer, ["Nxy", loads.get('Nxy', '')])
            w(writer, ["Mx", moments.get('Mx', '')])
            w(writer, ["My", moments.get('My', '')])
            w(writer, ["Mxy", moments.get('Mxy', '')])
            w(writer, ["DeltaT_degC", inputs.get('deltaT_degC', '')])

            # (Global allowables and derived strains removed; they are now per-lamina in the table above)

            w_blank(writer)
            # Optimization options
            w(writer, ["Optimization options"])
            opti = inputs.get('optimization', {})
            for k in [
                'enabled', 'SM_wanted', 'same_thickness_per_lam', 'total_thickness_target_mm',
                'lam_min_thick_mm', 'lam_max_thick_mm', 'criteria'
            ]:
                w(writer, [k, opti.get(k, '')])

            # Results section
            w_blank(writer, 2)
            w(writer, ["Results"])  # Section header
            if results is None:
                w(writer, ["(no results available)"])
            else:
                # Scalars and small vectors
                def write_named(name):
                    val = results.get(name)
                    if val is None:
                        return
                    if isinstance(val, list):
                        try:
                            # flatten only one level
                            w(writer, [name] + [str(x) for x in val])
                        except Exception:
                            w(writer, [name, str(val)])
                    else:
                        w(writer, [name, str(val)])

                # Engineering constants
                for key in [
                    'Ex', 'Ey', 'Gxy', 'nuxy', 'nuyx', 'Exb', 'Eyb', 'Dx', 'Dy', 't'
                ]:
                    write_named(key)

                # ABD matrices
                for mkey in ['A', 'B', 'D']:
                    mat = results.get(mkey)
                    if mat is not None:
                        w(writer, [mkey])
                        for row in mat:
                            w(writer, [" "] + list(map(str, row)))

                # z coordinates, eps0, kappa
                for key in ['z', 'eps0', 'kappa', 's1']:
                    write_named(key)

                # Per-ply lists
                def write_ply_list(name):
                    lst = results.get(name)
                    if lst is None:
                        return
                    w(writer, [name])
                    for i, item in enumerate(lst, start=1):
                        # item may be vector-like
                        if isinstance(item, (list, tuple)):
                            w(writer, [i] + [str(x) for x in item])
                        else:
                            w(writer, [i, str(item)])

                for key in [
                    'eps_i_global_list', 'stress_i_global_list',
                    'eps_i_local_list', 'stress_i_local_list',
                    'strain_FI_list', 'stress_FI_list'
                ]:
                    write_ply_list(key)

            # Assessment
            w_blank(writer)
            w(writer, ["Assessment"])
            if assess:
                w(writer, ["passed", assess.get('passed', '')])
                w(writer, ["min_sm_stress", assess.get('min_sm_stress', '')])
                w(writer, ["min_sm_strain", assess.get('min_sm_strain', '')])
                fails = assess.get('failures', [])
                if fails:
                    w(writer, ["failures (ply_index, type, component, FI)"])
                    for f in fails:
                        idx = f.get('index', '')
                        reasons = f.get('reasons', [])
                        for r in reasons:
                            # r = (type, comp, FI)
                            try:
                                rtype, comp, fi = r
                            except Exception:
                                rtype, comp, fi = ('', '', r)
                            w(writer, [idx, rtype, comp, fi])
                else:
                    w(writer, ["no failures"])
            else:
                w(writer, ["(no assessment)"])

            # Optimization result
            w_blank(writer)
            w(writer, ["Optimization result"])
            if optim:
                for k in ['success', 'message', 'nit', 'fun']:
                    if k in optim:
                        w(writer, [k, optim.get(k)])
                # thickness vectors
                if 'x' in optim:
                    w(writer, ['x_m'] + [str(v) for v in optim.get('x', [])])
                if 'x_mm' in optim:
                    w(writer, ['x_mm'] + [str(v) for v in optim.get('x_mm', [])])
                if 'min_SM_after' in optim:
                    w(writer, ['min_SM_after', optim.get('min_SM_after')])
                if 'total_thickness_after' in optim:
                    w(writer, ['total_thickness_after_m', optim.get('total_thickness_after')])
            else:
                w(writer, ["(no optimization result)"])


    def export_results_dialog(self):
        """Ask user where to save and export current results to CSV via external helper."""
        if not self.laminae_list:
            messagebox.showwarning("Export", "No laminae defined to export.")
            return
        # Ensure results exist
        if getattr(self, 'results', None) is None:
            try:
                self.update_results_display()
            except Exception:
                pass
        path = filedialog.asksaveasfilename(
            title="Save ICLC results",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if not path:
            return
        data = self.gather_export_data()
        try:
            self.export_iclc_results(path, data)
            messagebox.showinfo("Export", f"Results saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Export", f"Failed to export results.\n{e}")

### Method to run the application ###
    def run(self):
        """Start the main GUI loop"""
        self.root.mainloop()
    
# Start the application
print("\n" + "="*60)
print("Starting Interactive Composite Laminate Calculator...")
print("="*60)
app = LaminateCalculator() # create an object of the LaminateCalculator class
app.run()