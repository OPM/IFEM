// $Id$
//==============================================================================
//!
//! \file TestTimeStep.C
//!
//! \date Oct 8 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for the time stepping parameters.
//!
//==============================================================================

#include "TimeStep.h"

#include "Catch2Support.h"
#include "tinyxml2.h"

#include <vector>


namespace {

//! \brief Parses the time stepping parameters from an XML string.
void parse (TimeStep& tp, const char* xml)
{
  tinyxml2::XMLDocument doc;
  REQUIRE(doc.Parse(xml) == tinyxml2::XML_SUCCESS);
  REQUIRE(tp.parse(doc.RootElement()));
}

//! \brief Advances \a n steps and returns the times reached.
std::vector<double> advance (TimeStep& tp, int n)
{
  std::vector<double> times;
  for (int i = 0; i < n && tp.increment(); i++)
    times.push_back(tp.time.t);
  return times;
}

//! \brief Checks that the times reached are the expected ones.
void check (const std::vector<double>& times,
            const std::vector<double>& expected)
{
  REQUIRE(times.size() == expected.size());
  for (size_t i = 0; i < times.size(); i++)
    REQUIRE_THAT(times[i], WithinRel(expected[i]));
}

const char* twoSegments =
  "<timestepping maxCutbacks=\"3\">"
  "  <step start=\"0\" end=\"0.5\">0.25</step>"
  "  <step end=\"1.0\">0.1</step>"
  "</timestepping>";

}


TEST_CASE("TestTimeStep.Schedule")
{
  TimeStep tp;
  parse(tp,"<timestepping>"
           "  <step start=\"0\" end=\"1.0\">0.3</step>"
           "  <step end=\"2.0\">0.5</step>"
           "</timestepping>");
  REQUIRE(!tp.hasCutback());
  check(advance(tp,10),{
    0.3, 0.6, 0.9, 1.0, 1.5, 2.0});
  REQUIRE(tp.finished());
}


TEST_CASE("TestTimeStep.Cutback")
{
  TimeStep tp;
  parse(tp,"<timestepping maxCutbacks=\"2\">"
           "  <step start=\"0\" end=\"1.0\">0.25</step>"
           "</timestepping>");
  REQUIRE(tp.hasCutback());
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.t, WithinRel(0.0625));
  REQUIRE_THAT(tp.time.dt, WithinRel(0.0625));
  REQUIRE(tp.step == 1);

  // The step size grows back by f1 = 1.5 to the scheduled one
  check(advance(tp,5),{
    0.15625, 0.296875, 0.5078125, 0.7578125, 1.0});
  REQUIRE_THAT(tp.time.dt, WithinRel(0.2421875));
}


TEST_CASE("TestTimeStep.CutbackAtSegmentEnd")
{
  TimeStep tp;
  parse(tp,twoSegments);
  check(advance(tp,2),{
    0.25, 0.5});

  // The reduced step no longer reaches the end of the first segment,
  // which the next step is shortened to, and the second segment is intact
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.t, WithinRel(0.3125));
  check(advance(tp,4),{
    0.40625, 0.5, 0.6, 0.7});
}


TEST_CASE("TestTimeStep.MaxCutbacks")
{
  TimeStep tp;
  parse(tp,"<timestepping maxCutbacks=\"2\" f2=\"0.5\">"
           "  <step start=\"0\" end=\"1.0\">0.25</step>"
           "</timestepping>");
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.0625));
  REQUIRE(!tp.cutback());

  // The count is per step
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
}


TEST_CASE("TestTimeStep.CutbackMinStep")
{
  // Without a maximum number of cut-backs the minimum step size limits them
  TimeStep tp;
  parse(tp,"<timestepping dtMin=\"0.1\" dtMax=\"0.1\">"
           "  <step start=\"0\" end=\"1.0\">0.25</step>"
           "</timestepping>");
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.1));
  REQUIRE(!tp.cutback());

  // The following steps continue with the reduced step size
  check(advance(tp,2),{0.2, 0.3});

  TimeStep tp0;
  parse(tp0,"<timestepping>"
            "  <step start=\"0\" end=\"1.0\">0.25</step>"
            "</timestepping>");
  REQUIRE(tp0.increment());
  REQUIRE(!tp0.cutback());
}


TEST_CASE("TestTimeStep.CutbackConstantStep")
{
  TimeStep tp;
  parse(tp,"<timestepping start=\"0\" end=\"1.0\" dt=\"0.25\""
           " maxCutbacks=\"1\" f1=\"2\" f2=\"0.5\"/>");
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  check(advance(tp,4),{
    0.375, 0.625, 0.875, 1.0});
}


TEST_CASE("TestTimeStep.SuggestNextStep")
{
  TimeStep tp;
  parse(tp,twoSegments);
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.2);
  tp.suggestNextStep(0.1);
  tp.suggestNextStep(0.15);
  check(advance(tp,3),{
    0.35, 0.5, 0.6});

  // A suggestion larger than the scheduled step size has no effect
  tp.suggestNextStep(1.0);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.7));
}


TEST_CASE("TestTimeStep.SuggestConstantStep")
{
  // Without time step definitions, a suggestion replaces the step size
  TimeStep tp;
  parse(tp,"<timestepping start=\"0\" end=\"1.0\" dt=\"0.1\"/>");
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.3);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.4));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.7));
  tp.suggestNextStep(0.05);
  check(advance(tp,2),{0.75, 0.8});
}


TEST_CASE("TestTimeStep.SuggestMinStep")
{
  // A suggestion does not take the step size below the minimum one, with
  // and without time step definitions
  TimeStep tp;
  parse(tp,"<timestepping start=\"0\" end=\"1.0\" dt=\"0.2\""
           " dtMin=\"0.1\" dtMax=\"0.5\"/>");
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.01);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.3));

  TimeStep sched;
  parse(sched,"<timestepping dtMin=\"0.1\" dtMax=\"0.5\">"
              "  <step start=\"0\" end=\"1.0\">0.25</step>"
              "</timestepping>");
  REQUIRE(sched.increment());
  sched.suggestNextStep(0.01);
  REQUIRE(sched.increment());
  REQUIRE_THAT(sched.time.t, WithinRel(0.35));
}


TEST_CASE("TestTimeStep.SuggestCutback")
{
  // A cut-back reduces the step size to a smaller suggested one
  TimeStep tp;
  parse(tp,twoSegments);
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.2);
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.2));

  // but by the factor f2 if the suggestion is not smaller
  tp.suggestNextStep(0.3);
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.05));
}


TEST_CASE("TestTimeStep.TargetIterations")
{
  // The next step size is suggested from the most iterations reported in a
  // step, as the current one scaled by the target over them
  TimeStep tp;
  parse(tp,"<timestepping start=\"0\" end=\"1.0\" dt=\"0.1\""
           " target_iterations=\"3\" maxCutbacks=\"3\"/>");
  REQUIRE(tp.increment());
  tp.reportIterations(2);
  tp.reportIterations(6);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.05));

  // Iterations around those of the solvers override them
  tp.reportIterations(6);
  tp.setIterations(2);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.075));

  // Nothing is suggested without reported iterations, nor from those of a
  // step which is cut back
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.075));
  tp.reportIterations(1);
  REQUIRE(tp.cutback());
  const double dt = tp.time.dt;
  REQUIRE(tp.increment());
  REQUIRE(tp.time.dt <= dt*1.5);

  // A step trimmed to end a time step definition suggests from its size
  // before the trimming, so the next definition starts at its own size
  TimeStep trim;
  parse(trim,"<timestepping target_iterations=\"2\">"
             "  <step start=\"0\" end=\"0.25\">0.1</step>"
             "  <step start=\"0.25\" end=\"1\">0.2</step>"
             "</timestepping>");
  check(advance(trim,2),{0.1,0.2});
  trim.reportIterations(2);
  REQUIRE(trim.increment());
  REQUIRE_THAT(trim.time.dt, WithinRel(0.05));
  trim.reportIterations(1);
  REQUIRE(trim.increment());
  REQUIRE_THAT(trim.time.t, WithinRel(0.45));
  REQUIRE_THAT(trim.time.dt, WithinRel(0.2));

  // Without a target, the reported iterations are not used
  TimeStep plain;
  parse(plain,"<timestepping start=\"0\" end=\"1.0\" dt=\"0.1\"/>");
  REQUIRE(plain.increment());
  plain.reportIterations(1);
  REQUIRE(plain.increment());
  REQUIRE_THAT(plain.time.dt, WithinRel(0.1));
}


TEST_CASE("TestTimeStep.Copy")
{
  TimeStep tp;
  parse(tp,twoSegments);
  advance(tp,3);
  REQUIRE(tp.cutback());

  TimeStep copy(tp);
  check(advance(copy,4),advance(tp,4));
  REQUIRE(copy.step == tp.step);
}


#ifdef HAS_CEREAL
TEST_CASE("TestTimeStep.Serialize")
{
  TimeStep tp;
  parse(tp,twoSegments);
  advance(tp,2);
  REQUIRE(tp.cutback());
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.05);
  tp.iter = 12; // The iterations of the step, which adapt the next step size

  std::map<std::string,std::string> data;
  REQUIRE(tp.serialize(data));

  TimeStep restarted;
  parse(restarted,twoSegments);
  REQUIRE(restarted.deSerialize(data));
  REQUIRE(restarted.step == tp.step);
  REQUIRE(restarted.iter == 12);
  REQUIRE_THAT(restarted.time.t, WithinRel(tp.time.t));
  check(advance(restarted,5),advance(tp,5));

  // Restart data without the schedule state replays the time steps
  data.erase("TimeStep::schedule");
  TimeStep replayed;
  parse(replayed,twoSegments);
  REQUIRE(replayed.deSerialize(data));
  REQUIRE(replayed.step == 3);
  REQUIRE_THAT(replayed.time.t, WithinRel(0.6));
}


TEST_CASE("TestTimeStep.SerializeTargetIterations")
{
  // The iterations reported in the step, and the suggestion from them
  const char* target = "<timestepping start=\"0\" end=\"1.0\" dt=\"0.1\""
                       " target_iterations=\"3\"/>";
  TimeStep tp;
  parse(tp,target);
  REQUIRE(tp.increment());
  tp.reportIterations(6);
  std::map<std::string,std::string> data;
  REQUIRE(tp.serialize(data));
  TimeStep restarted;
  parse(restarted,target);
  REQUIRE(restarted.deSerialize(data));
  check(advance(restarted,3),advance(tp,3));
}


TEST_CASE("TestTimeStep.SerializeLongSegment")
{
  // The schedule state of a long step definition of constant size is
  // compact, and the restored time steps continue, also when replayed
  // from the start
  const char* xml = "<timestepping>"
                    "  <step start=\"0\" end=\"100\">0.01</step>"
                    "</timestepping>";
  TimeStep tp;
  parse(tp,xml);
  advance(tp,1000);

  // The size of the state does not grow with the steps
  std::map<std::string,std::string> data, later;
  REQUIRE(tp.serialize(data));
  TimeStep further(tp);
  advance(further,1000);
  REQUIRE(further.serialize(later));
  REQUIRE(later["TimeStep::schedule"].size() ==
          data["TimeStep::schedule"].size());

  TimeStep restarted;
  parse(restarted,xml);
  REQUIRE(restarted.deSerialize(data));
  check(advance(restarted,5),advance(tp,5));

  REQUIRE(restarted.reset(500));
  REQUIRE_THAT(restarted.time.t, WithinRel(5.0, 1.0e-12));
}


TEST_CASE("TestTimeStep.ReplayLimited")
{
  // A replay of the time steps reaches the times of the steps taken, also
  // where a suggestion or a cut-back limited them, and after a restart
  TimeStep tp;
  parse(tp,twoSegments);
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.1);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.35));
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.t, WithinRel(0.3875));
  REQUIRE(tp.increment());
  const double t4 = tp.time.t;

  std::map<std::string,std::string> data;
  REQUIRE(tp.serialize(data));
  TimeStep restarted;
  parse(restarted,twoSegments);
  REQUIRE(restarted.deSerialize(data));

  for (TimeStep* ts : { &tp, &restarted }) {
    REQUIRE(ts->reset(2));
    REQUIRE_THAT(ts->time.t, WithinRel(0.35));
    REQUIRE(ts->increment());
    REQUIRE_THAT(ts->time.t, WithinRel(0.3875));
    REQUIRE(ts->increment());
    REQUIRE_THAT(ts->time.t, WithinRel(t4));
  }
}


TEST_CASE("TestTimeStep.ReplayAdaptive")
{
  // The steps taken are replayed with their sizes, adapted to the
  // iterations when they were taken, while new steps are adapted
  TimeStep tp;
  parse(tp,"<timestepping dtMin=\"0.01\" dtMax=\"1.0\">"
           "  <step start=\"0\" end=\"10\">0.1</step>"
           "</timestepping>");
  const std::vector<double> taken = advance(tp,3);
  check(taken,{0.1, 0.2, 0.35});

  REQUIRE(tp.reset(3));
  REQUIRE_THAT(tp.time.t, WithinRel(0.35));
  REQUIRE(tp.reset(1));
  check(advance(tp,3),{0.2, 0.35, 0.575});

  // Also after a restart from a state behind the steps taken
  REQUIRE(tp.reset(1));
  std::map<std::string,std::string> data;
  REQUIRE(tp.serialize(data));
  TimeStep restarted;
  parse(restarted,"<timestepping dtMin=\"0.01\" dtMax=\"1.0\">"
                  "  <step start=\"0\" end=\"10\">0.1</step>"
                  "</timestepping>");
  REQUIRE(restarted.deSerialize(data));
  check(advance(restarted,3),{0.2, 0.35, 0.575});
}


TEST_CASE("TestTimeStep.RestartBehindSegment")
{
  // The recorded step sizes of a later definition are kept in a restart
  // from a state reset to an earlier one
  TimeStep tp;
  parse(tp,twoSegments);
  advance(tp,3);
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.t, WithinRel(0.525));
  REQUIRE(tp.reset(1));

  std::map<std::string,std::string> data;
  REQUIRE(tp.serialize(data));
  TimeStep restarted;
  parse(restarted,twoSegments);
  REQUIRE(restarted.deSerialize(data));
  check(advance(restarted,2),{0.5, 0.525});
}


TEST_CASE("TestTimeStep.ReplayCutRecovery")
{
  // A replay continues to grow the step size back after a cut-back, as the
  // steps taken did
  const char* xml = "<timestepping maxCutbacks=\"3\">"
                    "  <step start=\"0\" end=\"10\">0.25</step>"
                    "</timestepping>";
  TimeStep tp;
  parse(tp,xml);
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.09375));

  TimeStep uninterrupted(tp);
  REQUIRE(uninterrupted.increment());
  REQUIRE_THAT(uninterrupted.time.dt, WithinRel(0.140625));

  REQUIRE(tp.reset(2));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.140625));
  REQUIRE_THAT(tp.time.t, WithinRel(uninterrupted.time.t));
}


TEST_CASE("TestTimeStep.ResetConstantStep")
{
  // Without time step definitions, a reset to the current step keeps the
  // step size growing back after a cut-back
  const char* xml = "<timestepping start=\"0\" end=\"1\" dt=\"0.25\""
                    " maxCutbacks=\"3\"/>";
  TimeStep tp;
  parse(tp,xml);
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE(tp.increment());
  TimeStep uninterrupted(tp);
  REQUIRE(uninterrupted.increment());
  REQUIRE_THAT(uninterrupted.time.dt, WithinRel(0.140625));

  REQUIRE(tp.reset(2));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.140625));
}


TEST_CASE("TestTimeStep.ResetSuggestion")
{
  // A suggested size of the next step is kept in a reset to the current
  // step, and in a replay reaching it
  TimeStep tp;
  parse(tp,"<timestepping>"
           "  <step start=\"0\" end=\"10\">0.25</step>"
           "</timestepping>");
  REQUIRE(tp.increment());
  tp.suggestNextStep(0.1);
  REQUIRE(tp.reset(1));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.35));

  REQUIRE(tp.reset(1));
  tp.suggestNextStep(0.1);
  REQUIRE(tp.reset(0));
  check(advance(tp,2),{0.25, 0.35});
}


TEST_CASE("TestTimeStep.ReplayCutBack")
{
  // A cut-back of a replayed step recovers from it, as a new step would
  const char* xml = "<timestepping maxCutbacks=\"3\">"
                    "  <step start=\"0\" end=\"10\">0.25</step>"
                    "</timestepping>";
  TimeStep tp;
  parse(tp,xml);
  advance(tp,3);
  REQUIRE(tp.reset(1));
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.0625));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.09375));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.140625));
}


TEST_CASE("TestTimeStep.ReplayIterations")
{
  // A replay continues with the iterations of the steps taken, as the step
  // size is adapted to them
  TimeStep tp;
  parse(tp,"<timestepping dtMin=\"0.01\" dtMax=\"1.0\">"
           "  <step start=\"0\" end=\"10\">0.1</step>"
           "</timestepping>");
  for (int it : { 0, 12, 3 }) {
    REQUIRE(tp.increment());
    tp.iter = it;
  }
  TimeStep uninterrupted(tp);
  REQUIRE(uninterrupted.increment());

  REQUIRE(tp.reset(3));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(uninterrupted.time.t));
}


TEST_CASE("TestTimeStep.ReplayCutBackAbandons")
{
  // A cut-back of a replayed step abandons the steps taken after it
  TimeStep tp;
  parse(tp,"<timestepping maxCutbacks=\"3\">"
           "  <step start=\"0\" end=\"10\">0.25</step>"
           "</timestepping>");
  advance(tp,2);
  tp.suggestNextStep(0.01);
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.01));

  REQUIRE(tp.reset(1));
  REQUIRE(tp.increment());
  REQUIRE(tp.cutback());
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.09375));
}


TEST_CASE("TestTimeStep.ReplayCFL")
{
  // The step sizes of steps controlled by the CFL number are recorded too,
  // and a replay continues in the time step definitions as they did
  TimeStep tp;
  parse(tp,"<timestepping maxCFL=\"1.0\">"
           "  <step start=\"0\" end=\"1\">0.1 0.2</step>"
           "</timestepping>");
  REQUIRE(tp.increment());
  tp.time.CFL = 1.0;
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.t, WithinRel(0.25));
  tp.time.CFL = 0.0;
  TimeStep uninterrupted(tp);
  REQUIRE(uninterrupted.increment());
  REQUIRE_THAT(uninterrupted.time.dt, WithinRel(0.2));

  REQUIRE(tp.reset(1));
  REQUIRE_THAT(tp.time.t, WithinRel(0.1));
  REQUIRE(tp.reset(2));
  REQUIRE_THAT(tp.time.t, WithinRel(0.25));
  REQUIRE(tp.increment());
  REQUIRE_THAT(tp.time.dt, WithinRel(0.2));
  REQUIRE_THAT(tp.time.t, WithinRel(uninterrupted.time.t));
}


TEST_CASE("TestTimeStep.RestartLegacyAdaptive")
{
  // A restart without the schedule state replays the time step definitions
  // with the step sizes adapted, as when the steps were taken
  const char* xml = "<timestepping dtMin=\"0.01\" dtMax=\"1.0\">"
                    "  <step start=\"0\" end=\"10\">0.1</step>"
                    "</timestepping>";
  TimeStep tp;
  parse(tp,xml);
  advance(tp,3);
  REQUIRE_THAT(tp.time.t, WithinRel(0.35));

  std::map<std::string,std::string> data;
  REQUIRE(tp.serialize(data));
  data.erase("TimeStep::schedule");
  TimeStep restarted;
  parse(restarted,xml);
  REQUIRE(restarted.deSerialize(data));
  REQUIRE_THAT(restarted.time.t, WithinRel(0.35));
}
#endif
