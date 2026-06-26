#include "kinematics.hh"
#include <iostream>

template <int dim>
class CNCPathPhysics : public Kinematics<dim> {
private:
    double dtmin;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
private:
    struct LinearMove {
        Vector start;    // Absolute start position of the move
        Vector end;      // Absolute end position of the move
        double feed;     // Linear speed (same units as your system, e.g. mm/s)
        bool   rapid;    // True for rapid (G0) moves; currently treated same as feed moves
    };

    std::vector<LinearMove> m_moves;
    std::size_t             m_currentMove;  // index into m_moves

public:
    // --------------------------------------------------------------------
    // Constructor / destructor
    // --------------------------------------------------------------------
    CNCPathPhysics(NodeList* nodeList,
                   PhysicalConstants& constants)
        : Kinematics<dim>(nodeList, constants),
          dtmin(1e30),
          m_moves(),
          m_currentMove(0) {}

    virtual ~CNCPathPhysics() {}

    // --------------------------------------------------------------------
    // API for loading the toolpath from parsed G-code
    //
    // Typical usage from Python after parsing G0/G1:
    //   cncPhysics.addLinearMove(Vector({x, y, z}), feed, rapid);
    // --------------------------------------------------------------------
    void
    addLinearMove(const Vector& endPosition, double feed, bool rapid=false) {
        LinearMove move;
        move.feed  = feed;
        move.rapid = rapid;

        if (m_moves.empty()) {
            // First move: use current tool position from node 0 as the start
            NodeList* nodeList = this->nodeList;
            int numNodes       = nodeList->size();
            if (numNodes <= 0) {
                std::cerr << "CNCPathPhysics::addLinearMove: NodeList has no nodes!\n";
                // Start at origin as a fallback
                move.start = Vector();
            } else {
                VectorField* position = nodeList->template getField<Vector>("position");
                move.start            = position->getValue(0);
            }
        } else {
            // Subsequent moves start where the last one ended
            move.start = m_moves.back().end;
        }

        move.end = endPosition;
        m_moves.push_back(move);
    }

    void
    clearPath() {
        m_moves.clear();
        m_currentMove = 0;
    }

    bool
    pathComplete() const {
        return m_currentMove >= m_moves.size();
    }

    // --------------------------------------------------------------------
    // Physics interface overrides
    // --------------------------------------------------------------------

    // Name used in dt voting debug output
    virtual std::string
    name() const override {
        return "CNCPathPhysics";
    }

    // Estimate the next timestep based on how close we are to the end of the
    // current move. Integrator::VoteDt() will call this each cycle.
    virtual double
    EstimateTimestep() const override {
        return dtmin;
    }

    // Core derivative evaluator: for each node, set dx/dt so that the tool
    // moves along the current linear segment at the prescribed feed rate.
    virtual void
    EvaluateDerivatives(const State<dim>* initialState,
                        State<dim>&       deriv,
                        const double      time,
                        const double      dt) override {

        NodeList* nodeList = this->nodeList;
        int numNodes       = nodeList->size();

        VectorField* position = initialState->template getField<Vector>("position");
        VectorField* dxdt     = deriv       .template getField<Vector>("position");

        // We optionally support "velocity" as a diagnostic field on the State:
        VectorField* velocity = nullptr;
        VectorField* dvdt     = nullptr;
        try {
            velocity = initialState->template getField<Vector>("velocity");
            dvdt     = deriv       .template getField<Vector>("velocity");
        } catch (...) {
            // If no velocity field is present, we just won't set dvdt
        }

        const double bigNumber = 1e30;
        double       local_dtmin = bigNumber;

        // No path loaded or already finished: zero derivatives, huge dt
        if (m_moves.empty() || m_currentMove >= m_moves.size()) {
            if (dxdt)  dxdt->fill(1,Vector());
            if (dvdt)  dvdt->fill(1,Vector());
            dtmin = bigNumber;
            this->lastDt = dt;
            return;
        }

        const LinearMove& move = m_moves[m_currentMove];

        // Basic direction of the whole move (start -> end)
        Vector fullDir = move.end - move.start;
        double fullLen2 = fullDir.mag2();
        double fullLen  = (fullLen2 > 0.0 ? std::sqrt(fullLen2) : 0.0);

        // Guard against degenerate move
        if (fullLen == 0.0 || move.feed <= 0.0) {
            if (dxdt)  dxdt->fill(1,Vector());
            if (dvdt)  dvdt->fill(1,Vector());
            dtmin = bigNumber;
            this->lastDt = dt;
            return;
        }

        // Unit direction of the move
        Vector unitDir = fullDir * (1.0 / fullLen);

        // For CNC, we typically have one node, but we'll loop for generality
        for (int i = 0; i < numNodes; ++i) {
            Vector x = position->getValue(i);

            // Remaining distance to end of this move
            Vector diff = move.end - x;
            double remaining2 = diff.mag2();
            double remaining  = (remaining2 > 0.0 ? std::sqrt(remaining2) : 0.0);

            // If we've effectively reached the end of this move, advance the
            // move index; on the *next* call we'll pick up the new segment.
            const double eps = 1e-9;
            if (remaining < eps) {
                // Snap position exactly to end on this node in the NodeList
                // (not strictly required, but helps keep the path tight)
                VectorField* nodePos = nodeList->template getField<Vector>("position");
                nodePos->setValue(i, move.end);

                // Only advance move index once (for node 0) to avoid racing ahead
                if (i == 0 && m_currentMove + 1 < m_moves.size()) {
                    ++m_currentMove;
                }

                // No motion this substep
                dxdt->setValue(i, Vector());
                if (dvdt) dvdt->setValue(i, Vector());
                continue;
            }

            // Desired speed along the path
            double speed = move.feed;

            // Direction we actually move along this step is toward the end point
            Vector stepDir = diff * (1.0 / remaining);  // normalized diff

            // dx/dt = v = speed * direction
            Vector v = stepDir * speed;
            dxdt->setValue(i, v);

            // If we keep velocity as a state field, then dv/dt = 0 for kinematic motion
            if (dvdt) {
                dvdt->setValue(i, Vector());  // no acceleration in this simple model
            }

            // Time to reach the end of this move from this node
            double tToEnd = remaining / speed;
            local_dtmin = std::min(local_dtmin, tToEnd);
        }

        // If for some reason local_dtmin never changed, keep dt huge to avoid
        // constraining other packages.
        if (local_dtmin == bigNumber) {
            dtmin = bigNumber;
        } else {
            dtmin = local_dtmin;
        }

        this->lastDt = dt;
    }
};