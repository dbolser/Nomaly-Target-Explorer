-- Create users table (if not exists)
DROP TABLE IF EXISTS users2;
CREATE TABLE IF NOT EXISTS users2 (
    id INT AUTO_INCREMENT PRIMARY KEY,
    username VARCHAR(100) UNIQUE NOT NULL,
    password VARCHAR(255) NOT NULL,
    email VARCHAR(100),
    is_active BOOLEAN DEFAULT TRUE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Create simplified user_permissions table
DROP TABLE IF EXISTS user_permissions;
CREATE TABLE IF NOT EXISTS user_permissions (
    id INT AUTO_INCREMENT PRIMARY KEY,
    user_id INT NOT NULL,
    allowed_paths TEXT NOT NULL,  -- Comma-separated list of allowed phecodes, or '*' for admin
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (user_id) REFERENCES users2(id),
    UNIQUE KEY unique_user (user_id)  -- One permission entry per user
);

