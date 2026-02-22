#version 330 core
//phong lighting
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;

uniform vec3 lightPos;
uniform vec3 viewPos;

void main() {
    vec3 norm= normalize(Normal);
    vec3 lightColor = vec3(1.0);
    vec3 objectColor = vec3(0.3,0.6,1.0); // Light blue

    //Ambient 
    float ambientStrength = 0.2;
    vec3 ambient = ambientStrength * lightColor ;

    //Diffuse
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    //Specular 
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos -FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength*spec*lightColor;

    // ===== Dark Rim (輪郭を少し暗くする) =====
    float rim = 1.0 - max(dot(norm, viewDir), 0.0);
    rim = pow(rim, 1.5);   // 強さ調整

    // 暗くする係数（0.0〜0.3くらいで調整）
    float darken = 0.25 * rim;

    vec3 result = (ambient + diffuse + specular)*objectColor;
    // ★ 縁だけ少し暗く
    result *= (1.0 - darken);
    FragColor = vec4(result,1.0);
}
