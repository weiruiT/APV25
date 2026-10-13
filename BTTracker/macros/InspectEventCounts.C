#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>

namespace {

void InspectOneTree(TTree *tree, const char *cut)
{
    const Long64_t treeEntries = tree->GetEntries();

    std::cout << "\nTree: " << tree->GetName() << "\n";
    std::cout << "  tree entries       : " << treeEntries << "\n";

    if (!tree->GetBranch("event")) {
        std::cout << "  event branch       : not present\n";
        return;
    }

    // 保证TTree::Draw保存全部结果，而不只保存默认Estimate数量。
    const Long64_t oldEstimate = tree->GetEstimate();
    tree->SetEstimate(treeEntries + 1);

    const Long64_t selected =
        tree->Draw("event", cut ? cut : "", "goff");

    if (selected < 0) {
        std::cout << "  ERROR              : invalid cut\n";
        tree->SetEstimate(oldEstimate);
        return;
    }

    std::set<Long64_t> eventIDs;
    const double *values = tree->GetV1();

    for (Long64_t i = 0; i < selected; ++i) {
        eventIDs.insert(
            static_cast<Long64_t>(std::llround(values[i]))
        );
    }

    tree->SetEstimate(oldEstimate);

    std::cout << "  selection          : "
              << ((cut && cut[0]) ? cut : "(none)") << "\n";
    std::cout << "  selected entries   : " << selected << "\n";
    std::cout << "  unique events      : " << eventIDs.size() << "\n";

    if (eventIDs.empty())
        return;

    const Long64_t minEvent = *eventIDs.begin();
    const Long64_t maxEvent = *eventIDs.rbegin();
    const Long64_t span = maxEvent - minEvent + 1;
    const Long64_t missing =
        span - static_cast<Long64_t>(eventIDs.size());
    const Long64_t repeated =
        selected - static_cast<Long64_t>(eventIDs.size());

    std::cout << "  event ID range     : ["
              << minEvent << ", " << maxEvent << "]\n";
    std::cout << "  event ID span      : " << span << "\n";
    std::cout << "  missing IDs in span: " << missing << "\n";
    std::cout << "  repeated entries   : " << repeated << "\n";
    std::cout << "  entries/event      : "
              << std::fixed << std::setprecision(3)
              << static_cast<double>(selected) / eventIDs.size()
              << "\n";
}

} // namespace

void InspectEventCounts(const char *fileName,
                        const char *treeName = "",
                        const char *cut = "")
{
    TFile file(fileName, "READ");

    if (file.IsZombie()) {
        std::cerr << "Cannot open file: " << fileName << "\n";
        return;
    }

    std::cout << "File: " << file.GetName() << "\n";

    // 指定tree：可以使用cut。
    if (treeName && treeName[0]) {
        TTree *tree = nullptr;
        file.GetObject(treeName, tree);

        if (!tree) {
            std::cerr << "Tree not found: " << treeName << "\n";
            return;
        }

        InspectOneTree(tree, cut);
        return;
    }

    // 未指定tree：检查文件内所有顶层tree。
    if (cut && cut[0]) {
        std::cerr
            << "A tree name must be specified when using a cut.\n";
        return;
    }

    std::set<std::string> inspectedNames;
    TIter next(file.GetListOfKeys());

    while (TKey *key = static_cast<TKey *>(next())) {
        const std::string name = key->GetName();

        // 同一个tree可能有多个cycle，只检查最新cycle。
        if (!inspectedNames.insert(name).second)
            continue;

        TObject *object = file.Get(name.c_str());
        TTree *tree = dynamic_cast<TTree *>(object);

        if (tree)
            InspectOneTree(tree, "");
    }
}
