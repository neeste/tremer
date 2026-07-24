import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import subprocess
import threading
import os
import sys
import webbrowser

# Add the 'tools' directory to the path so we can import update_tremer_ages
sys.path.append(os.path.join(os.path.dirname(__file__), 'tools'))
import update_tremer_ages

class RedirectText:
    """A helper class to redirect stdout/stderr to a Tkinter Text widget."""
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END) # Auto-scroll to the bottom

    def flush(self):
        pass

class TremerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("TREMER Control Panel")
        self.root.geometry("600x500")
        
        self.project_file = None
        self.output_files = []

        # --- UI LAYOUT ---
        top_frame = tk.Frame(root, pady=10)
        top_frame.pack(fill=tk.X, padx=10)

        # File Selection
        self.file_label = tk.Label(top_frame, text="No project selected", font=("Helvetica", 12, "italic"), fg="gray")
        self.file_label.pack(side=tk.LEFT, expand=True, fill=tk.X)

        btn_select = tk.Button(top_frame, text="Select Project File...", command=self.select_file, font=("Helvetica", 12))
        btn_select.pack(side=tk.RIGHT)

        # Options Frame
        opt_frame = tk.Frame(root, pady=5)
        opt_frame.pack(fill=tk.X, padx=10)
        
        self.relaxed_var = tk.BooleanVar()
        chk_relaxed = tk.Checkbutton(opt_frame, text="Heuristic Discovery Mode (Relaxed - Allow distance > 1)", variable=self.relaxed_var)
        chk_relaxed.pack(side=tk.LEFT)

        # Actions Frame
        action_frame = tk.Frame(root, pady=10)
        action_frame.pack(fill=tk.X, padx=10)

        self.btn_update = tk.Button(action_frame, text="1. Update FTDNA Ages", command=self.run_update, state=tk.DISABLED, font=("Helvetica", 12))
        self.btn_update.pack(side=tk.LEFT, padx=5)

        self.btn_run = tk.Button(action_frame, text="2. Run TREMER Engine", command=self.run_engine, state=tk.DISABLED, font=("Helvetica", 12))
        self.btn_run.pack(side=tk.LEFT, padx=5)
        
        self.btn_view = tk.Button(action_frame, text="3. View Output HTML", command=self.view_output, state=tk.DISABLED, font=("Helvetica", 12))
        self.btn_view.pack(side=tk.RIGHT, padx=5)

        # Log Output Console
        log_frame = tk.Frame(root)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        tk.Label(log_frame, text="Console Output:", anchor="w").pack(fill=tk.X)
        self.log_text = scrolledtext.ScrolledText(log_frame, bg="black", fg="lightgreen", font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # Redirect standard output to the text widget
        sys.stdout = RedirectText(self.log_text)
        sys.stderr = RedirectText(self.log_text)
        print("Welcome to TREMER. Please select a project text file to begin.")

    def select_file(self):
        filepath = filedialog.askopenfilename(
            title="Select TREMER Project File",
            filetypes=(("Text Files", "*.txt"), ("All Files", "*.*"))
        )
        if filepath:
            self.project_file = filepath
            filename = os.path.basename(filepath)
            self.file_label.config(text=f"Selected: {filename}", fg="black", font=("Helvetica", 12, "bold"))
            self.btn_update.config(state=tk.NORMAL)
            self.btn_run.config(state=tk.NORMAL)
            print(f"\n--- Loaded Project: {filepath} ---")

    def run_update(self):
        if not self.project_file: return
        self.btn_update.config(state=tk.DISABLED)
        self.btn_run.config(state=tk.DISABLED)
        
        def task():
            try:
                print("\n[Updating FTDNA Ages...]")
                update_tremer_ages.process_tremer_file(self.project_file)
                print("[Update Complete!]")
            except Exception as e:
                print(f"[ERROR] {e}")
            finally:
                self.btn_update.config(state=tk.NORMAL)
                self.btn_run.config(state=tk.NORMAL)
        
        threading.Thread(target=task, daemon=True).start()

    def run_engine(self):
        if not self.project_file: return
        self.btn_update.config(state=tk.DISABLED)
        self.btn_run.config(state=tk.DISABLED)
        self.btn_view.config(state=tk.DISABLED)
        
        def task():
            try:
                print("\n[Running TREMER Engine...]")
                
                # Make sure the engine exists, if not run make build
                engine_path = os.path.join(os.path.dirname(__file__), "tremer_rewrite")
                if not os.path.exists(engine_path):
                    print("Compiling engine...")
                    subprocess.run(["make", "build"], cwd=os.path.dirname(__file__), check=True)
                
                cmd = [engine_path, self.project_file]
                if self.relaxed_var.get():
                    cmd.append("-relaxed")
                
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=os.path.dirname(__file__))
                for line in process.stdout:
                    print(line, end="")
                process.wait()
                
                print(f"[Engine Complete with exit code {process.returncode}]")
                if process.returncode == 0:
                    self.detect_outputs()
            except Exception as e:
                print(f"[ERROR] {e}")
            finally:
                self.btn_update.config(state=tk.NORMAL)
                self.btn_run.config(state=tk.NORMAL)
        
        threading.Thread(target=task, daemon=True).start()

    def detect_outputs(self):
        import glob
        base_dir = os.path.dirname(__file__)
        html_files = glob.glob(os.path.join(base_dir, "*_b*.html")) + glob.glob(os.path.join(base_dir, "*_g*.html"))
        
        # Sort by modification time to find the newest ones
        html_files.sort(key=os.path.getmtime, reverse=True)
        
        if html_files:
            self.output_files = html_files
            self.btn_view.config(state=tk.NORMAL)
            basename = os.path.basename(html_files[0])
            self.btn_view.config(text=f"3. View Output ({basename})")
            print(f"Detected outputs: {basename}")

    def view_output(self):
        if self.output_files:
            filepath = os.path.abspath(self.output_files[0])
            print(f"Opening browser for: {filepath}")
            webbrowser.open(f"file://{filepath}")

if __name__ == "__main__":
    root = tk.Tk()
    app = TremerApp(root)
    root.mainloop()
