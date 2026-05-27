---
description: Work on a GitHub ticket end-to-end in a Worktree
allowed-tools: Read, Write, Edit, Glob, Grep, Bash(git:*), Bash(gh:*), Bash(npm:*), mcp__github__*
---

# Ticket Workflow

Work on ticket: $ARGUMENTS

## Instructions

### 1. Read the Ticket

First, fetch and understand the ticket:

```
Use the GitHub MCP tool to:
- Get ticket details (title, description, acceptance criteria)
- Check linked tickets
- Review any comments or attachments
```

Summarize:
- What needs to be done
- Acceptance criteria
- Any blockers or dependencies

### 2. Explore the Codebase

Before coding:
- Search for related code
- Understand the current implementation
- Identify files that need changes

### 3. Create a Branch

We always branch from staging for development.
The branch name should be a combination of username,
ticket ID and a brief description: -

    {initials}-{ticket-id}-{brief-description}

Then create a Worktree using the repository name for the branch: -

```bash
git fetch origin
git worktree add -b {repository-name}_{branch} ../{branch} staging
cd ../{repository-name}_{branch}
```

### 4. Implement the Changes

- Follow project patterns (check relevant skills)
- Write tests first (TDD)
- Make incremental commits

### 5. Update the Ticket

As you work:
- Add comments with progress updates
- Update status (In Progress → In Review)
- Log any blockers or questions

### 6. Create PR and Link

When ready:
- Create PR with `gh pr create`
- Link the PR to the ticket
- Add ticket ID to PR title: `feat(123): description`

### 7. If You Find a Bug

If you discover an unrelated bug while working:
1. Create a new ticket with details
2. Link it to the current ticket if related
3. Note it in the PR description
4. Continue with original task

## Example Workflow

```
Me: /ticket 123

Claude:
1. Fetching 123 from GitHub...
   Title: Add user profile avatar upload
   Description: Users should be able to upload a profile picture...
   Acceptance Criteria:
   - [ ] Upload button on profile page
   - [ ] Support JPG/PNG up to 5MB
   - [ ] Show loading state during upload

2. Searching codebase for profile-related code...
   Found: src/screens/Profile/ProfileScreen.tsx
   Found: src/components/Avatar/Avatar.tsx

3. Creating worktree: fragalysis-backend_cw-123-avatar-upload

4. [Implements feature with TDD approach]

5. Creating PR and linking to 123...
   PR #456 created: feat(123): add avatar upload to profile
```
