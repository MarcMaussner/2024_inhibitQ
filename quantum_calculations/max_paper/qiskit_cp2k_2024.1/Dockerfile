FROM cp2k/cp2k:2024.1_openmpi_generic_psmp

RUN apt update && apt install -y python3 python3-pip git
RUN pip3 install --upgrade pip
RUN pip3 install qiskit qiskit-aer qiskit-algorithms

RUN git clone -b fix-qiskit1.x-compatibility https://github.com/KarimElgammal/qiskit-nature-cp2k.git && \
    cd qiskit-nature-cp2k && \
    pip3 install .

WORKDIR /workspace

RUN echo '#!/bin/bash\ntail -f /dev/null' > /keep_alive.sh && chmod +x /keep_alive.sh
ENTRYPOINT ["/keep_alive.sh"]